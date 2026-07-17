# =============================================================================
# RBMRB -- ChemicalShiftStatistics module
#
# CONFIRMED search API response format (PyBMRB 3.0.9 source):
#   GET /v2/search/chemical_shifts?comp_id=ALA&atom_id=CA
#   Returns: {
#     "columns": ["Atom_chem_shift.Val", "Sample_conditions.pH",
#                 "Sample_conditions.Temperature_K", "Atom_chem_shift.Comp_ID",
#                 "Atom_chem_shift.Atom_ID", "Atom_chem_shift.Ambiguity_code",
#                 "Atom_chem_shift.Entry_ID", "Atom_chem_shift.Entity_assembly_ID",
#                 "Atom_chem_shift.Comp_index_ID", ...],
#     "data": [["52.1", "7.0", "298", "ALA", "CA", "1", "15060", "1", "5"], ...]
#   }
# =============================================================================

.STANDARD_AAS <- c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
                    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP",
                    "TYR","VAL")
.STANDARD_NAS <- c("A","C","G","U","DA","DC","DG","DT")

#' Fetch chemical shift data from the BMRB database
#'
#' Downloads chemical shifts for the given residue and atom combination from
#' the BMRB database. Mirrors `pybmrb.ChemicalShiftStatistics.get_data()`.
#'
#' @param residue Residue in IUPAC 3-letter format (e.g. "ALA"). "*" = any.
#' @param atom    Atom name (e.g. "CA", "HB*"). "*" = any.
#' @param filtered  Remove outliers beyond sd_limit * SD. Default TRUE.
#' @param sd_limit  SD factor for outlier removal. Default 10.
#' @param ambiguity Ambiguity code filter ("*" = no filter; 1 = unambiguous).
#' @param ph_min,ph_max  pH range filter (optional).
#' @param t_min,t_max    Temperature range in Kelvin (optional).
#' @param standard_amino_acids Restrict to standard amino acids. Default TRUE.
#' @param database "macromolecules" (default) or "metabolomics".
#' @return A data.frame with columns Val, Comp_ID, Atom_ID, Ambiguity_code,
#'   Entry_ID, Comp_index_ID, pH, Temperature_K and others.
#' @export
#' @examples
#' \dontrun{
#'   df <- get_cs_data("ALA", "CA")
#'   df <- get_cs_data("*", "HB*")
#' }
get_cs_data <- function(residue = "*",
                         atom    = "*",
                         filtered = TRUE,
                         sd_limit = 10,
                         ambiguity = "*",
                         ph_min = NULL, ph_max = NULL,
                         t_min  = NULL, t_max  = NULL,
                         standard_amino_acids = TRUE,
                         database = "macromolecules") {

  database <- match.arg(database, c("macromolecules","metabolomics"))

  # Build params -- only comp_id and atom_id are sent to the API
  # pH and temperature are filtered client-side (matching PyBMRB behaviour)
  params <- list()
  if (!is.null(residue) && residue != "*") params$comp_id  <- residue
  if (!is.null(atom)    && atom    != "*") params$atom_id  <- atom

  res_label <- if (is.null(residue) || residue == "*") "all residues" else residue
  atm_label <- if (is.null(atom)    || atom    == "*") "all atoms"    else atom
  cli::cli_progress_message("Querying BMRB for {res_label} / {atm_label} ...")

  raw <- tryCatch(
    .bmrb_get("/search/chemical_shifts", params = params),
    error = function(e) {
      cli::cli_alert_danger("API call failed: {conditionMessage(e)}")
      NULL
    }
  )
  if (is.null(raw)) return(data.frame())

  # The API may return either:
  #   {"columns": [...], "data": [...]}  (newer format)
  #   {"tags":    [...], "data": [...]}  (older format)
  # Both use dotted names like "Atom_chem_shift.Val" for columns/tags.
  col_key <- names(raw)[names(raw) %in% c("columns", "tags")]
  if (length(col_key) == 0L || is.null(raw$data)) {
    cli::cli_alert_warning(
      "Unexpected response keys: {paste(names(raw), collapse=',')}.")
    return(data.frame())
  }
  col_key <- col_key[[1L]]

  cols <- unlist(raw[[col_key]], use.names = FALSE)
  rows <- raw$data

  if (length(rows) == 0L) {
    cli::cli_alert_warning("BMRB returned 0 rows for {res_label}/{atm_label}.")
    return(data.frame())
  }

  mat <- do.call(rbind, lapply(rows, function(r)
    vapply(r,
           function(v) if (is.null(v)||identical(v,".")||identical(v,"?"))
                         NA_character_ else as.character(v),
           character(1L))))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- cols

  # Shorten column names: "Atom_chem_shift.Val" -> "Val" etc.
  colnames(df) <- sub("^Atom_chem_shift\\.", "", colnames(df))
  colnames(df) <- sub("^Sample_conditions\\.", "", colnames(df))
  colnames(df) <- sub("^Assigned_chem_shift_list\\.", "", colnames(df))

  # Numeric coercions
  for (col in c("Val","Val_err","pH","Temperature_K"))
    if (col %in% colnames(df))
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  if ("Ambiguity_code" %in% colnames(df))
    df[["Ambiguity_code"]] <- suppressWarnings(as.integer(df[["Ambiguity_code"]]))
  if ("Comp_index_ID" %in% colnames(df))
    df[["Comp_index_ID"]] <- suppressWarnings(as.integer(df[["Comp_index_ID"]]))

  # Client-side filters (matching PyBMRB exactly)
  if (standard_amino_acids && (is.null(residue) || residue == "*") &&
      "Comp_ID" %in% colnames(df))
    df <- df[df$Comp_ID %in% c(.STANDARD_AAS, .STANDARD_NAS), ]

  if (!identical(ambiguity, "*") && "Ambiguity_code" %in% colnames(df))
    df <- df[!is.na(df$Ambiguity_code) &
               df$Ambiguity_code == as.integer(ambiguity), ]

  ph_col <- if ("pH" %in% colnames(df)) "pH" else NULL
  if (!is.null(ph_col) && !is.null(ph_min) && !is.null(ph_max))
    df <- df[!is.na(df[[ph_col]]) &
               df[[ph_col]] >= ph_min & df[[ph_col]] <= ph_max, ]

  t_col <- if ("Temperature_K" %in% colnames(df)) "Temperature_K" else NULL
  if (!is.null(t_col) && !is.null(t_min) && !is.null(t_max))
    df <- df[!is.na(df[[t_col]]) &
               df[[t_col]] >= t_min & df[[t_col]] <= t_max, ]

  if (filtered && "Val" %in% colnames(df) && nrow(df) > 0L)
    df <- .filter_outliers_df(df, "Val", sd_limit)

  cli::cli_alert_info("{nrow(df)} shifts loaded for {res_label}/{atm_label}.")
  df
}

#' Fetch paired 2D chemical shift data for two atoms in the same residue
#'
#' Mirrors `pybmrb.ChemicalShiftStatistics.get_2d_chemical_shifts()`.
#'
#' @param residue  Residue name (IUPAC 3-letter).
#' @param atom1    First atom (X axis).
#' @param atom2    Second atom (Y axis).
#' @param filtered,sd_limit  See [get_cs_data()].
#' @param ambiguity1,ambiguity2  Per-atom ambiguity filters.
#' @param ph_min,ph_max,t_min,t_max  Condition filters.
#' @return Named list with `atom1` and `atom2` numeric vectors.
#' @export
#' @examples \dontrun{ paired <- get_2d_cs_data("ALA", "CA", "CB") }
get_2d_cs_data <- function(residue, atom1, atom2,
                            filtered = TRUE, sd_limit = 10,
                            ambiguity1 = "*", ambiguity2 = "*",
                            ph_min = NULL, ph_max = NULL,
                            t_min  = NULL, t_max  = NULL) {
  df1 <- get_cs_data(residue, atom1, filtered=filtered, sd_limit=sd_limit,
                      ambiguity=ambiguity1,
                      ph_min=ph_min, ph_max=ph_max, t_min=t_min, t_max=t_max)
  df2 <- get_cs_data(residue, atom2, filtered=filtered, sd_limit=sd_limit,
                      ambiguity=ambiguity2,
                      ph_min=ph_min, ph_max=ph_max, t_min=t_min, t_max=t_max)

  if (nrow(df1) == 0L || nrow(df2) == 0L) {
    cli::cli_alert_warning("No data for one or both atoms.")
    return(list(atom1 = numeric(0), atom2 = numeric(0)))
  }

  key_cols <- intersect(c("Entry_ID","Comp_index_ID","Entity_assembly_ID"),
                         intersect(colnames(df1), colnames(df2)))
  if (length(key_cols) == 0L) {
    cli::cli_alert_warning("Cannot merge: no shared key columns.")
    return(list(atom1 = numeric(0), atom2 = numeric(0)))
  }

  merged <- merge(df1[, c(key_cols,"Val")],
                  df2[, c(key_cols,"Val")],
                  by = key_cols, suffixes = c("_1","_2"))
  merged <- merged[!is.na(merged$Val_1) & !is.na(merged$Val_2), ]
  list(atom1 = merged$Val_1, atom2 = merged$Val_2)
}

#' Fetch raw BMRB chemical shift data (no outlier filtering)
#'
#' Mirrors `pybmrb.ChemicalShiftStatistics.get_data_from_bmrb()`.
#'
#' @inheritParams get_cs_data
#' @param list_of_atoms Optional "RES-ATOM" vector, e.g. `c("ALA-CB","GLY-CA")`.
#' @return A data.frame of raw chemical shift records.
#' @export
#' @examples
#' \dontrun{
#'   df  <- get_cs_from_bmrb_db("ALA", "CB")
#'   df2 <- get_cs_from_bmrb_db(list_of_atoms = c("ALA-CB","GLY-CA"))
#' }
get_cs_from_bmrb_db <- function(residue = "*", atom = "*",
                                 list_of_atoms = NULL,
                                 ambiguity = "*",
                                 ph_min = NULL, ph_max = NULL,
                                 t_min = NULL, t_max = NULL,
                                 standard_amino_acids = TRUE,
                                 database = "macromolecules") {
  if (!is.null(list_of_atoms)) {
    dfs <- lapply(list_of_atoms, function(pair) {
      parts <- strsplit(pair, "-", fixed = TRUE)[[1L]]
      if (length(parts) != 2L)
        cli::cli_abort("list_of_atoms must be 'RESIDUE-ATOM', got: {pair}")
      get_cs_from_bmrb_db(parts[[1L]], parts[[2L]],
                           ambiguity=ambiguity,
                           ph_min=ph_min, ph_max=ph_max,
                           t_min=t_min, t_max=t_max,
                           database=database)
    })
    dfs <- Filter(Negate(is.null), dfs)
    if (length(dfs) == 0L) return(data.frame())
    return(do.call(rbind, dfs))
  }
  get_cs_data(residue=residue, atom=atom, filtered=FALSE,
               ambiguity=ambiguity,
               ph_min=ph_min, ph_max=ph_max,
               t_min=t_min, t_max=t_max,
               standard_amino_acids=standard_amino_acids,
               database=database)
}

# Internal
.filter_outliers_df <- function(df, col, sd_limit) {
  vals <- df[[col]]
  mu   <- mean(vals, na.rm = TRUE)
  sig  <- stats::sd(vals, na.rm = TRUE)
  if (is.na(sig) || sig == 0) return(df)
  keep <- !is.na(vals) & vals >= mu - sd_limit*sig & vals <= mu + sd_limit*sig
  df[keep, , drop = FALSE]
}

.standardise_cs_cols <- function(cols) {
  map <- c(
    entry_id              = "Entry_ID",
    entity_assembly_id    = "Entity_assembly_ID",
    comp_index_id         = "Comp_index_ID",
    comp_id               = "Comp_ID",
    atom_id               = "Atom_ID",
    atom_type             = "Atom_type",
    val                   = "Val",
    val_err               = "Val_err",
    ambiguity_code        = "Ambiguity_code",
    assigned_chem_shift_list_id = "Assigned_chem_shift_list_ID",
    ph                    = "pH",
    temperature           = "Temperature",
    temperature_k         = "Temperature_K"
  )
  lc <- tolower(cols)
  for (i in seq_along(cols))
    if (lc[[i]] %in% names(map)) cols[[i]] <- map[[lc[[i]]]]
  cols
}
