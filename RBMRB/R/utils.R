# =============================================================================
# RBMRB -- shared utility functions (no duplicate .output_figure here)
# =============================================================================

#' List the 20 standard amino acid codes
#' @return Character vector of 3-letter amino acid codes.
#' @export
#' @examples standard_amino_acids()
standard_amino_acids <- function() .STANDARD_AAS

#' List standard nucleic acid residue codes
#' @return Character vector of RNA and DNA residue codes.
#' @export
#' @examples standard_nucleic_acids()
standard_nucleic_acids <- function() .STANDARD_NAS

#' Summarise a nested chemical shift list
#'
#' Prints a tidy summary of datasets, chains, and atom coverage.
#'
#' @param cs_data Nested list from [cs_from_bmrb()] or [cs_from_file()].
#' @return Invisibly, a data.frame summary.
#' @export
#' @examples
#' \dontrun{ cs <- cs_from_bmrb(15060); cs_summary(cs) }
cs_summary <- function(cs_data) {
  rows <- list()
  for (ds in names(cs_data)) {
    chains <- setdiff(names(cs_data[[ds]]), "seq_ids")
    for (ch in chains) {
      cd      <- cs_data[[ds]][[ch]]
      seq_ids <- cd$seq_ids
      atoms   <- character(0)
      for (sid in as.character(seq_ids))
        atoms <- union(atoms, names(cd[[sid]]))
      rows <- c(rows, list(data.frame(
        Dataset    = ds, Chain = ch,
        N_residues = length(seq_ids),
        Seq_start  = if (length(seq_ids)) min(seq_ids) else NA_integer_,
        Seq_end    = if (length(seq_ids)) max(seq_ids) else NA_integer_,
        Atoms      = paste(sort(atoms), collapse=","),
        stringsAsFactors = FALSE)))
    }
  }
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  print(df)
  invisible(df)
}

#' Convert a nested CS list to a tidy data.frame
#'
#' @param cs_data Nested list from [cs_from_bmrb()] or [cs_from_file()].
#' @return A data.frame with columns Dataset, Chain, Seq_ID, Atom_ID, Val.
#' @export
#' @examples
#' \dontrun{ cs <- cs_from_bmrb(15060); df <- cs_to_df(cs) }
cs_to_df <- function(cs_data) {
  rows <- list()
  for (ds in names(cs_data)) {
    chains <- setdiff(names(cs_data[[ds]]), "seq_ids")
    for (ch in chains) {
      cd      <- cs_data[[ds]][[ch]]
      seq_ids <- cd$seq_ids
      for (sid in seq_ids) {
        rd <- cd[[as.character(sid)]]
        if (is.null(rd)) next
        for (atm in names(rd))
          rows <- c(rows, list(data.frame(
            Dataset=ds, Chain=ch, Seq_ID=sid, Atom_ID=atm,
            Val=rd[[atm]], stringsAsFactors=FALSE)))
      }
    }
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

#' Remove outliers from a numeric vector
#'
#' @param x        Numeric vector.
#' @param sd_limit Number of standard deviations beyond which to remove.
#' @return Numeric vector with outliers removed.
#' @export
#' @examples filter_cs_outliers(c(52, 55, 56, 54, 300), sd_limit = 2)
filter_cs_outliers <- function(x, sd_limit = 10) {
  x   <- x[!is.na(x)]
  if (length(x) < 2L) return(x)
  mu  <- mean(x)
  sig <- stats::sd(x)
  if (is.na(sig) || sig == 0) return(x)
  x[x >= mu - sd_limit * sig & x <= mu + sd_limit * sig]
}

#' Compute summary statistics for a set of chemical shifts
#'
#' @param x Numeric vector of chemical shift values.
#' @return A named numeric vector: n, mean, sd, median, min, max.
#' @export
#' @examples cs_stats(c(55.1, 55.3, 56.0, 54.8))
cs_stats <- function(x) {
  x <- x[!is.na(x)]
  c(n      = length(x),
    mean   = mean(x),
    sd     = stats::sd(x),
    median = stats::median(x),
    min    = min(x),
    max    = max(x))
}

#' Peek at a nested CS structure
#'
#' @param cs_data  Output of [cs_from_bmrb()] or [cs_from_file()].
#' @param n_res    Number of residues to show per chain.
#' @export
#' @examples
#' \dontrun{ cs <- cs_from_bmrb(15060); cs_peek(cs) }
cs_peek <- function(cs_data, n_res = 5L) {
  for (ds in names(cs_data)) {
    cat(sprintf("\nDataset: %s\n", ds))
    for (ch in setdiff(names(cs_data[[ds]]), "seq_ids")) {
      cd   <- cs_data[[ds]][[ch]]
      sids <- utils::head(cd$seq_ids, n_res)
      cat(sprintf("  Chain %s  (first %d residues)\n", ch, length(sids)))
      for (sid in sids) {
        atoms <- cd[[as.character(sid)]]
        cat(sprintf("    Seq %d: %s\n", sid,
          paste(names(atoms), round(unlist(atoms), 2L),
                sep="=", collapse="  ")))
      }
    }
  }
  invisible(cs_data)
}
