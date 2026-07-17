# =============================================================================
# RBMRB -- NMR-STAR 3 file parser
# Parses the Atom_chem_shift loop from .str / .nmrstar files.
# This is a pure-R parser; no external dependencies needed.
# =============================================================================

#' Parse one or more NMR-STAR files and extract Atom_chem_shift data
#'
#' @param file_paths Character vector of NMR-STAR file paths.
#' @param auth_tag   If `TRUE`, use `Auth_seq_ID` instead of `Comp_index_ID`.
#' @param data_set_ids Optional character vector of dataset labels (one per file).
#' @return A list of data.frames, one per file, with chemical shift columns.
#' @noRd
.parse_nmrstar_files <- function(file_paths, auth_tag = FALSE,
                                  data_set_ids = NULL) {
  out <- vector("list", length(file_paths))
  for (i in seq_along(file_paths)) {
    fp  <- file_paths[[i]]
    did <- if (!is.null(data_set_ids)) data_set_ids[[i]] else basename(fp)
    if (!file.exists(fp)) {
      cli::cli_alert_warning("File not found, skipping: {fp}")
      next
    }
    cli::cli_progress_message("Parsing {fp} ...")
    raw <- readLines(fp, warn = FALSE, encoding = "UTF-8")
    dfs <- .extract_cs_loops(raw)
    if (length(dfs) == 0L) {
      cli::cli_alert_warning("No Atom_chem_shift loop found in {fp}")
      next
    }
    df <- do.call(rbind, dfs)
    df[["_dataset_id"]] <- did
    out[[i]] <- df
  }
  out[!vapply(out, is.null, logical(1L))]
}

# ---------------------------------------------------------------------------
# Core loop extractor
# ---------------------------------------------------------------------------

.extract_cs_loops <- function(lines) {
  n    <- length(lines)
  loops_out <- list()
  i    <- 1L

  while (i <= n) {
    line <- trimws(lines[[i]])

    # ---- Detect start of a loop_ block ------------------------------------
    if (tolower(line) == "loop_") {
      i <- i + 1L
      # Collect tag names
      tags <- character(0)
      while (i <= n) {
        tline <- trimws(lines[[i]])
        if (startsWith(tline, "_")) {
          tags <- c(tags, tline)
          i    <- i + 1L
        } else {
          break
        }
      }

      # Is this an Atom_chem_shift loop?
      is_cs_loop <- any(startsWith(tags, "_Atom_chem_shift."))
      short_tags <- sub("^_Atom_chem_shift\\.", "", tags)

      # Read data rows until stop_
      rows <- list()
      while (i <= n) {
        dline <- trimws(lines[[i]])
        if (tolower(dline) == "stop_") {
          i <- i + 1L
          break
        }
        # Skip blank and comment lines
        if (nchar(dline) == 0L || startsWith(dline, "#")) {
          i <- i + 1L
          next
        }
        # Tokenise the line (handles quoted strings)
        tokens <- .tokenise_nmrstar(dline)
        if (length(tokens) > 0L && is_cs_loop) {
          rows <- c(rows, list(tokens))
        }
        i <- i + 1L
      }

      if (is_cs_loop && length(rows) > 0L) {
        # Some rows may span multiple lines; accumulate until we have ntags tokens
        ntags  <- length(short_tags)
        merged <- list()
        buf    <- character(0)
        for (r in rows) {
          buf <- c(buf, r)
          while (length(buf) >= ntags) {
            merged <- c(merged, list(buf[seq_len(ntags)]))
            buf    <- buf[-seq_len(ntags)]
          }
        }
        # Build data.frame
        if (length(merged) > 0L) {
          mat <- do.call(rbind, merged)
          df  <- as.data.frame(mat, stringsAsFactors = FALSE)
          colnames(df) <- short_tags[seq_len(ncol(df))]
          # Type coercions
          df <- .coerce_cs_df(df)
          loops_out <- c(loops_out, list(df))
        }
      }
    } else {
      i <- i + 1L
    }
  }
  loops_out
}

# ---------------------------------------------------------------------------
# Tokeniser: splits a line respecting single-quoted strings
# ---------------------------------------------------------------------------

.tokenise_nmrstar <- function(line) {
  # Handles: plain tokens, 'quoted strings', and . (null)
  tokens <- character(0)
  i      <- 1L
  nc     <- nchar(line)
  while (i <= nc) {
    ch <- substr(line, i, i)
    if (ch == " " || ch == "\t") {
      i <- i + 1L
    } else if (ch == "'") {
      # Find closing quote not followed by non-whitespace
      j <- i + 1L
      repeat {
        if (j > nc) break
        if (substr(line, j, j) == "'") {
          if (j == nc || substr(line, j + 1L, j + 1L) %in% c(" ", "\t")) break
        }
        j <- j + 1L
      }
      tokens <- c(tokens, substr(line, i + 1L, j - 1L))
      i <- j + 1L
    } else if (ch == "\"") {
      j <- i + 1L
      while (j <= nc && substr(line, j, j) != "\"") j <- j + 1L
      tokens <- c(tokens, substr(line, i + 1L, j - 1L))
      i <- j + 1L
    } else {
      j <- i + 1L
      while (j <= nc && !substr(line, j, j) %in% c(" ", "\t")) j <- j + 1L
      tokens <- c(tokens, substr(line, i, j - 1L))
      i <- j
    }
  }
  tokens
}

# ---------------------------------------------------------------------------
# Type coercions for Atom_chem_shift columns
# ---------------------------------------------------------------------------

.coerce_cs_df <- function(df) {
  na_vals <- c(".", "?", "NA", "")
  for (col in c("Val", "Val_err", "Assign_fig_of_merit")) {
    if (col %in% colnames(df)) {
      df[[col]] <- suppressWarnings(as.numeric(ifelse(df[[col]] %in% na_vals, NA, df[[col]])))
    }
  }
  for (col in c("Comp_index_ID", "Seq_ID", "Ambiguity_code",
                "Entity_assembly_ID", "Entity_ID")) {
    if (col %in% colnames(df)) {
      df[[col]] <- suppressWarnings(as.integer(ifelse(df[[col]] %in% na_vals, NA, df[[col]])))
    }
  }
  # Replace remaining "." with NA in character columns
  char_cols <- setdiff(colnames(df), c("Val", "Val_err", "Assign_fig_of_merit",
                                        "Comp_index_ID", "Seq_ID",
                                        "Ambiguity_code", "Entity_assembly_ID",
                                        "Entity_ID"))
  for (col in char_cols) {
    df[[col]] <- ifelse(df[[col]] %in% na_vals, NA_character_, df[[col]])
  }
  df
}

# ---------------------------------------------------------------------------
# Determine the chain column to use
# ---------------------------------------------------------------------------

.chain_col <- function(df, auth_tag = FALSE) {
  if (auth_tag && "Auth_asym_ID" %in% colnames(df) &&
      !all(is.na(df[["Auth_asym_ID"]]))) {
    return("Auth_asym_ID")
  }
  if ("Entity_assembly_ID" %in% colnames(df)) return("Entity_assembly_ID")
  NA_character_
}

.seq_col <- function(df, auth_tag = FALSE) {
  if (auth_tag && "Auth_seq_ID" %in% colnames(df) &&
      !all(is.na(df[["Auth_seq_ID"]]))) {
    return("Auth_seq_ID")
  }
  if ("Comp_index_ID" %in% colnames(df)) return("Comp_index_ID")
  "Seq_ID"
}
