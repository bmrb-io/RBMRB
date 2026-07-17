# =============================================================================
# RBMRB 1.1.0 — ChemicalShift module
# Uses .fetch_entry_json() which fetches ?format=zlib (PyBMRB method)
# =============================================================================

#' Extract chemical shifts from BMRB database entries
#'
#' Fetches assigned chemical shift data for one or more BMRB entries.
#' Uses the same zlib-compressed JSON endpoint as PyBMRB / pynmrstar.
#' Mirrors `pybmrb.ChemicalShift.from_bmrb()`.
#'
#' @param bmrb_ids A single BMRB entry ID or a vector of IDs.
#' @param auth_tag Logical. Use `Auth_seq_ID` for sequence numbering. Default FALSE.
#' @return A nested named list:
#'   `list(dataset_id = list(chain = list(seq_id = list(atom = cs_value),
#'   seq_ids = integer_vector)))`.
#' @export
#' @examples
#' \dontrun{
#'   cs <- cs_from_bmrb(15060)
#'   cs_peek(cs)
#'   cs_from_bmrb(c(17076, 17077))
#' }
cs_from_bmrb <- function(bmrb_ids, auth_tag = FALSE) {
  ids    <- as.character(.ensure_vector(bmrb_ids))
  result <- list()

  for (id in ids) {
    cli::cli_progress_message("Fetching BMRB entry {id} ...")

    entry_json <- tryCatch(
      .fetch_entry_json(id),
      error = function(e) {
        cli::cli_alert_warning("Could not fetch entry {id}: {conditionMessage(e)}")
        NULL
      }
    )
    if (is.null(entry_json)) next

    dfs <- .entry_json_to_cs_dfs(entry_json)
    if (length(dfs) == 0L) {
      cli::cli_alert_warning(
        "No chemical shift data found for entry {id}. ",
        "Run bmrb_debug_entry({id}) for diagnostics.")
      next
    }

    df      <- do.call(rbind, dfs)
    n_valid <- sum(!is.na(df$Val))
    cli::cli_alert_info("Entry {id}: {n_valid} chemical shifts loaded.")

    nested <- .cs_df_to_nested(df, dataset_id = id, auth_tag = auth_tag)
    result <- c(result, nested)
  }
  result
}

#' Extract chemical shifts from local NMR-STAR files
#'
#' Mirrors `pybmrb.ChemicalShift.from_file()`.
#'
#' @param file_paths Character vector of NMR-STAR file paths.
#' @param auth_tag   Logical. Use `Auth_seq_ID`. Default FALSE.
#' @param data_set_ids Optional character vector of dataset labels.
#' @return Same nested-list structure as [cs_from_bmrb()].
#' @export
#' @examples
#' \dontrun{
#'   cs <- cs_from_file("my_protein.str")
#' }
cs_from_file <- function(file_paths, auth_tag = FALSE, data_set_ids = NULL) {
  file_paths <- as.character(.ensure_vector(file_paths))
  if (!is.null(data_set_ids)) {
    data_set_ids <- as.character(data_set_ids)
    if (length(data_set_ids) != length(file_paths))
      cli::cli_abort("`data_set_ids` must have the same length as `file_paths`.")
  }
  raw_dfs <- .parse_nmrstar_files(file_paths, auth_tag = auth_tag,
                                   data_set_ids = data_set_ids)
  result <- list()
  for (df in raw_dfs) {
    did    <- unique(df[["_dataset_id"]])[[1L]]
    df_sub <- df[, setdiff(colnames(df), "_dataset_id"), drop = FALSE]
    nested <- .cs_df_to_nested(df_sub, dataset_id = did, auth_tag = auth_tag)
    result <- c(result, nested)
  }
  result
}

#' Get a flat data.frame of chemical shifts for one or more BMRB entries
#'
#' @inheritParams cs_from_bmrb
#' @return A data.frame suitable for use with ggplot2 or dplyr.
#' @export
#' @examples
#' \dontrun{ df <- cs_table(c(15060, 17076)) }
cs_table <- function(bmrb_ids, auth_tag = FALSE) {
  ids <- as.character(.ensure_vector(bmrb_ids))
  dfs <- lapply(ids, function(id) {
    ej <- tryCatch(
      .fetch_entry_json(id),
      error = function(e) {
        cli::cli_alert_warning("Failed {id}: {conditionMessage(e)}"); NULL })
    if (is.null(ej)) return(NULL)
    raw_dfs <- .entry_json_to_cs_dfs(ej)
    if (length(raw_dfs) == 0L) return(NULL)
    df <- do.call(rbind, raw_dfs)
    df[["Entry_ID"]] <- id
    df
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0L) return(data.frame())
  do.call(rbind, dfs)
}
