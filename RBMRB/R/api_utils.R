# =============================================================================
# RBMRB -- BMRB API v2 helpers and public utility functions
# =============================================================================

.BMRB_API <- "https://api.bmrb.io/v2"
.RBMRB_UA <- "RBMRB/3.0.0 (https://github.com/uwbmrb/RBMRB)"

# ---------------------------------------------------------------------------
# Internal HTTP helpers
# ---------------------------------------------------------------------------

.bmrb_get_raw <- function(path, params = list(), timeout = 90) {
  url <- paste0(.BMRB_API, path)
  if (length(params) > 0L) {
    qs <- paste(mapply(function(k, v)
      paste0(utils::URLencode(as.character(k), reserved = TRUE), "=",
             utils::URLencode(as.character(v), reserved = TRUE)),
      names(params), params, USE.NAMES = FALSE), collapse = "&")
    url <- paste0(url, "?", qs)
  }
  req  <- httr2::request(url) |>
    httr2::req_headers(Application = .RBMRB_UA) |>
    httr2::req_timeout(timeout)
  resp <- tryCatch(httr2::req_perform(req),
    error = function(e) cli::cli_abort("BMRB API error: {conditionMessage(e)}"))
  if (httr2::resp_status(resp) != 200L)
    cli::cli_abort("HTTP {httr2::resp_status(resp)} from BMRB ({path}).")
  httr2::resp_body_raw(resp)
}

.bmrb_get <- function(path, params = list(), timeout = 60) {
  jsonlite::fromJSON(rawToChar(.bmrb_get_raw(path, params, timeout)),
                     simplifyVector = FALSE)
}

.fetch_entry_json <- function(bmrb_id) {
  id        <- as.character(bmrb_id)
  raw_bytes <- tryCatch(.bmrb_get_raw(paste0("/entry/", id),
                                       params = list(format = "zlib")),
                         error = function(e) NULL)
  if (!is.null(raw_bytes) && length(raw_bytes) > 0L) {
    for (dtype in c("unknown", "gzip")) {
      txt <- tryCatch(rawToChar(memDecompress(raw_bytes, type = dtype)),
                      error = function(e) NULL)
      if (!is.null(txt) && nchar(txt) > 10L) {
        ej <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                       error = function(e) NULL)
        if (!is.null(ej) && !is.null(ej[["saveframes"]])) return(ej)
        break
      }
    }
  }
  ej <- tryCatch(.bmrb_get(paste0("/entry/", id)), error = function(e) NULL)
  if (is.null(ej)) return(NULL)
  if (!"saveframes" %in% names(ej)) {
    ekey <- setdiff(names(ej), c("error","error_message","error_type","bmrb_id"))
    if (length(ekey) >= 1L) {
      inner <- ej[[ekey[[1L]]]]
      if (is.list(inner)) {
        if ("saveframes" %in% names(inner)) {
          ej <- inner
        } else {
          converted <- lapply(inner, function(sf) {
            lps <- list()
            for (nm in names(sf))
              for (lp in sf[[nm]]) {
                lp[["category"]] <- paste0("_", nm)
                lps <- c(lps, list(lp))
              }
            list(loops = lps)
          })
          ej <- list(saveframes = converted, bmrb_id = id)
        }
      }
    }
  }
  ej
}

.loop_to_df <- function(loop_json) {
  tags <- unlist(loop_json$tags, use.names = FALSE)
  rows <- loop_json$data
  if (length(rows) == 0L) {
    df <- as.data.frame(matrix(character(0L), 0L, length(tags)),
                        stringsAsFactors = FALSE)
    colnames(df) <- tags; return(df)
  }
  mat <- do.call(rbind, lapply(rows, function(r)
    vapply(r, function(v)
      if (is.null(v) || identical(v, ".") || identical(v, "?"))
        NA_character_ else as.character(v), character(1L))))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- tags
  for (col in c("Val", "Val_err"))
    if (col %in% colnames(df))
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  if ("Ambiguity_code" %in% colnames(df))
    df[["Ambiguity_code"]] <- suppressWarnings(as.integer(df[["Ambiguity_code"]]))
  for (col in c("Comp_index_ID", "Seq_ID", "Entity_assembly_ID"))
    if (col %in% colnames(df))
      df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  df
}

.entry_json_to_cs_dfs <- function(entry_json) {
  out        <- list()
  saveframes <- entry_json[["saveframes"]]
  if (is.null(saveframes) || length(saveframes) == 0L) return(out)
  for (sf in saveframes) {
    loops <- sf[["loops"]]
    if (is.null(loops) || length(loops) == 0L) next
    for (lp in loops) {
      cat_name <- lp[["category"]]
      if (!is.null(cat_name) &&
          grepl("Atom_chem_shift", cat_name, ignore.case = TRUE) &&
          !is.null(lp$tags) && !is.null(lp$data) && length(lp$data) > 0L) {
        df <- .loop_to_df(lp)
        if (nrow(df) > 0L) out <- c(out, list(df))
      }
    }
  }
  out
}

.cs_df_to_nested <- function(df, dataset_id, auth_tag = FALSE) {
  seq_col <- if (auth_tag && "Auth_seq_ID" %in% colnames(df) &&
                   !all(is.na(df[["Auth_seq_ID"]])))
               "Auth_seq_ID" else "Comp_index_ID"
  chain_col <- if ("Auth_asym_ID" %in% colnames(df) &&
                     !all(is.na(df[["Auth_asym_ID"]])))
                 "Auth_asym_ID" else "Entity_assembly_ID"
  needed <- c(seq_col, "Atom_ID", "Val")
  if (!all(needed %in% colnames(df))) return(list())
  if (!chain_col %in% colnames(df)) { df[["._ch"]] <- "A"; chain_col <- "._ch" }
  df[[chain_col]] <- ifelse(is.na(df[[chain_col]]) | df[[chain_col]] == ".",
                             "A", df[[chain_col]])
  result <- list()
  for (chain in unique(df[[chain_col]])) {
    cdf        <- df[df[[chain_col]] == chain & !is.na(df$Val), ]
    seq_ids    <- sort(unique(cdf[[seq_col]]))
    chain_data <- list(seq_ids = seq_ids)
    for (sid in seq_ids) {
      sdf <- cdf[!is.na(cdf[[seq_col]]) & cdf[[seq_col]] == sid, ]
      chain_data[[as.character(sid)]] <- as.list(setNames(sdf$Val, sdf$Atom_ID))
    }
    result[[as.character(chain)]] <- chain_data
  }
  setNames(list(result), as.character(dataset_id))
}

.null_na       <- function(x, d = NA_character_)
  if (is.null(x) || length(x) == 0L) d else as.character(x[[1L]])
.ensure_vector <- function(x)
  if (length(x) == 1L && is.list(x)) unlist(x) else x
.is_metabolomics <- function(id)
  grepl("^bmse", as.character(id), ignore.case = TRUE)

# ---------------------------------------------------------------------------
# Public API utility functions
# ---------------------------------------------------------------------------

#' Print the installed RBMRB version
#'
#' Prints the package version and lists all exported functions.
#'
#' @return The version string, invisibly.
#' @export
#' @examples
#' RBMRB_version()
RBMRB_version <- function() {
  cat("RBMRB version:", as.character(utils::packageVersion("RBMRB")), "\n")
  invisible(TRUE)
}

#' Run a connectivity self-test
#'
#' Performs three checks: BMRB API status, entry fetch and parse for entry
#' 15060, and the chemical shift search endpoint. Useful for verifying
#' that the BMRB API is reachable from your network.
#'
#' @return A named logical vector with elements \code{api_status},
#'   \code{entry_fetch}, and \code{search_query} (invisibly).
#' @export
#' @examples
#' \dontrun{
#' bmrb_test()
#' }
bmrb_test <- function() {
  cat("RBMRB", as.character(utils::packageVersion("RBMRB")), "self-test\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  results <- logical(3)

  cat("[1/3] API status   ... ")
  results[1] <- tryCatch(!is.null(.bmrb_get("/status")), error = function(e) FALSE)
  cat(if (results[1]) "PASS\n" else "FAIL\n")

  cat("[2/3] Entry 15060  ... ")
  results[2] <- tryCatch({
    ej  <- .fetch_entry_json("15060")
    dfs <- .entry_json_to_cs_dfs(ej)
    n   <- sum(vapply(dfs, nrow, integer(1L)))
    if (n > 0L) { cat("PASS (", n, " shifts)\n", sep = ""); TRUE }
    else { cat("FAIL\n"); FALSE }
  }, error = function(e) { cat("FAIL\n"); FALSE })

  cat("[3/3] Search ALA CA ... ")
  results[3] <- tryCatch({
    r <- .bmrb_get("/search/chemical_shifts",
                    params = list(comp_id = "ALA", atom_id = "CA"))
    k <- names(r)[names(r) %in% c("columns", "tags")]
    if (length(k) > 0L && length(r$data) > 0L) {
      cat("PASS (", length(r$data), " rows)\n", sep = ""); TRUE
    } else { cat("FAIL\n"); FALSE }
  }, error = function(e) { cat("FAIL\n"); FALSE })

  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(if (all(results)) "All tests PASSED.\n" else "Some tests FAILED.\n")
  names(results) <- c("api_status", "entry_fetch", "search_query")
  invisible(results)
}

#' Diagnose BMRB entry fetch problems
#'
#' Fetches a BMRB entry using both the zlib-compressed and plain JSON
#' endpoints and reports what is found at each step. Run this when
#' \code{\link{cs_from_bmrb}} returns no data.
#'
#' @param bmrb_id A single BMRB entry ID (integer or character).
#' @return Invisibly, the parsed entry JSON list, or \code{NULL} on failure.
#' @export
#' @examples
#' \dontrun{
#' bmrb_debug_entry(15060)
#' }
bmrb_debug_entry <- function(bmrb_id) {
  id <- as.character(bmrb_id)
  cat("=== bmrb_debug_entry(", id, ") ===\n")
  cat("[1] Fetching zlib bytes ... ")
  raw <- tryCatch(
    .bmrb_get_raw(paste0("/entry/", id), params = list(format = "zlib")),
    error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); NULL })
  if (!is.null(raw)) {
    cat(length(raw), "bytes\n")
    for (dt in c("unknown", "gzip")) {
      txt <- tryCatch(rawToChar(memDecompress(raw, type = dt)),
                      error = function(e) NULL)
      if (!is.null(txt) && nchar(txt) > 10L) {
        cat("[2] Decompressed (type =", dt, "):", nchar(txt), "chars\n")
        ej <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                       error = function(e) NULL)
        if (!is.null(ej) && !is.null(ej$saveframes)) {
          dfs <- .entry_json_to_cs_dfs(ej)
          n   <- sum(vapply(dfs, nrow, integer(1L)))
          cat("[3] CS rows:", n, "--", if (n > 0L) "SUCCESS\n" else "FAIL\n")
          return(invisible(ej))
        }
        break
      }
    }
    cat("[2] Decompression failed\n")
  }
  cat("[1b] Plain JSON fallback ... ")
  ej <- tryCatch(.bmrb_get(paste0("/entry/", id)),
                  error = function(e) { cat("FAIL\n"); NULL })
  if (!is.null(ej)) {
    cat("keys:", paste(names(ej), collapse = ", "), "\n")
    ej <- .fetch_entry_json(id)
    if (!is.null(ej)) {
      dfs <- .entry_json_to_cs_dfs(ej)
      n   <- sum(vapply(dfs, nrow, integer(1L)))
      cat("[3] CS rows:", n, "--", if (n > 0L) "SUCCESS\n" else "FAIL\n")
      return(invisible(ej))
    }
  }
  cat("FAIL: could not fetch entry", id, "\n")
  invisible(NULL)
}

#' Get BMRB API server status
#'
#' Returns current server information including entry counts and
#' database last-update timestamps from the BMRB REST API.
#'
#' @return A list with server status information.
#' @export
#' @examples
#' \dontrun{
#' bmrb_status()
#' }
bmrb_status <- function() .bmrb_get("/status")

#' List all BMRB entry IDs
#'
#' Returns all public entry IDs in the specified BMRB database.
#'
#' @param database One of \code{"macromolecules"} (default),
#'   \code{"metabolomics"}, or \code{"chemcomps"}.
#' @return A character vector of BMRB entry IDs.
#' @export
#' @examples
#' \dontrun{
#' ids <- bmrb_list_entries()
#' length(ids)
#' }
bmrb_list_entries <- function(database = "macromolecules") {
  database <- match.arg(database, c("macromolecules", "metabolomics", "chemcomps"))
  unlist(.bmrb_get("/list_entries", params = list(database = database)),
         use.names = FALSE)
}

#' Search BMRB using the instant-search endpoint
#'
#' Searches BMRB by title, author name, organism, amino-acid sequence,
#' InChI string, or other terms.
#'
#' @param term A search string.
#' @param database One of \code{"macromolecules"}, \code{"metabolomics"},
#'   or \code{"combined"} (default).
#' @return A data.frame with columns \code{value}, \code{label},
#'   \code{sub_date}, \code{authors}, and \code{link}.
#' @export
#' @examples
#' \dontrun{
#' bmrb_search("ubiquitin human")
#' bmrb_search("GB3", database = "macromolecules")
#' }
bmrb_search <- function(term, database = "combined") {
  res <- .bmrb_get("/instant", params = list(term = term, database = database))
  if (length(res) == 0L) return(data.frame())
  do.call(rbind, lapply(res, function(e)
    data.frame(value    = .null_na(e$value),
               label    = .null_na(e$label),
               sub_date = .null_na(e$sub_date),
               authors  = paste(unlist(e$authors), collapse = "; "),
               link     = .null_na(e$link),
               stringsAsFactors = FALSE)))
}

#' Fetch raw metadata for a BMRB entry
#'
#' @param bmrb_ids A single BMRB entry ID or a character/integer vector.
#' @return A named list of raw entry metadata (API JSON as R list).
#' @export
#' @examples
#' \dontrun{
#' meta <- bmrb_fetch_entry(15060)
#' names(meta[["15060"]])
#' }
bmrb_fetch_entry <- function(bmrb_ids) {
  ids <- as.character(.ensure_vector(bmrb_ids))
  setNames(lapply(ids, function(id) .fetch_entry_json(id)), ids)
}

#' Get the citation for a BMRB entry
#'
#' @param bmrb_id A single BMRB entry ID.
#' @param format Output format: \code{"bibtex"} (default), \code{"json-ld"},
#'   or \code{"text"}.
#' @return A character string containing the citation.
#' @export
#' @examples
#' \dontrun{
#' cat(bmrb_citation(15060, format = "text"))
#' }
bmrb_citation <- function(bmrb_id, format = "bibtex") {
  format <- match.arg(format, c("bibtex", "json-ld", "text"))
  raw <- .bmrb_get(paste0("/entry/", bmrb_id, "/citation"),
                    params = list(format = format))
  if (is.character(raw)) raw else paste(unlist(raw), collapse = "\n")
}

#' Get PDB IDs associated with a BMRB entry
#'
#' @param bmrb_id A single BMRB entry ID.
#' @return A data.frame with columns \code{pdb_id}, \code{match_type},
#'   and \code{comment}.
#' @export
#' @examples
#' \dontrun{
#' bmrb_pdb_ids(15000)
#' }
bmrb_pdb_ids <- function(bmrb_id) {
  res <- .bmrb_get(paste0("/search/get_pdb_ids_from_bmrb_id/", bmrb_id))
  if (length(res) == 0L) return(data.frame())
  do.call(rbind, lapply(res, function(r)
    data.frame(pdb_id     = .null_na(r$pdb_id),
               match_type = .null_na(r$match_type),
               comment    = .null_na(r$comment),
               stringsAsFactors = FALSE)))
}

#' Get BMRB IDs associated with a PDB entry
#'
#' @param pdb_id A PDB ID string (e.g. \code{"2JM0"}).
#' @return A data.frame with columns \code{bmrb_id}, \code{match_type},
#'   and \code{comment}.
#' @export
#' @examples
#' \dontrun{
#' bmrb_from_pdb("2JM0")
#' }
bmrb_from_pdb <- function(pdb_id) {
  res <- .bmrb_get(paste0("/search/get_bmrb_ids_from_pdb_id/", pdb_id))
  if (length(res) == 0L) return(data.frame())
  do.call(rbind, lapply(res, function(r)
    data.frame(bmrb_id    = .null_na(r$bmrb_id),
               match_type = .null_na(r$match_type),
               comment    = .null_na(r$comment),
               stringsAsFactors = FALSE)))
}

#' Upload an NMR-STAR file to BMRB for temporary access
#'
#' Uploads a local NMR-STAR file to the BMRB server. The entry is
#' stored for 7 days. The returned token can be used as a BMRB ID
#' in all other package functions.
#'
#' @param file_path Path to a local NMR-STAR \code{.str} file.
#' @return A character string token usable as a BMRB ID (invisibly).
#' @export
#' @examples
#' \dontrun{
#' token <- bmrb_upload_entry("my_protein.str")
#' n15hsqc(bmrb_ids = token)
#' }
bmrb_upload_entry <- function(file_path) {
  if (!file.exists(file_path)) cli::cli_abort("File not found: {file_path}")
  body <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  req  <- httr2::request(paste0(.BMRB_API, "/entry/")) |>
    httr2::req_headers(Application    = .RBMRB_UA,
                       `Content-Type` = "text/plain") |>
    httr2::req_method("POST") |>
    httr2::req_body_raw(chartr("\r", "", body), type = "text/plain")
  res  <- httr2::resp_body_json(httr2::req_perform(req), simplifyVector = FALSE)
  cli::cli_alert_info("Uploaded. Token: {res[['entry_id']]}")
  invisible(res[["entry_id"]])
}
