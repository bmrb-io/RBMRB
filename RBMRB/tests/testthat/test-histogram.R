# =============================================================================
# tests/testthat/test-histogram.R
# Tests for histogram module and chemical shift statistics
# =============================================================================

# ---------------------------------------------------------------------------
# ChemicalShiftStatistics internal helpers
# ---------------------------------------------------------------------------

test_that(".filter_outliers_df removes values beyond sd_limit", {
  df      <- data.frame(Val = c(52.1, 52.5, 51.9, 53.0, 52.3, 300.0),
                        stringsAsFactors = FALSE)
  filtered <- RBMRB:::.filter_outliers_df(df, "Val", sd_limit = 2)
  expect_s3_class(filtered, "data.frame")
  expect_true(nrow(filtered) < nrow(df))
  expect_false(any(filtered$Val == 300.0))
})

test_that(".standardise_cs_cols maps lowercase to standard names", {
  cols <- c("comp_id", "atom_id", "val", "entry_id", "temperature")
  std  <- RBMRB:::.standardise_cs_cols(cols)
  expect_equal(std, c("Comp_ID","Atom_ID","Val","Entry_ID","Temperature"))
})

test_that(".STANDARD_AAS has exactly 20 entries", {
  expect_length(RBMRB:::.STANDARD_AAS, 20L)
})

test_that(".STANDARD_NAS covers RNA and DNA", {
  nas <- RBMRB:::.STANDARD_NAS
  expect_true("A"  %in% nas)  # RNA adenine
  expect_true("DT" %in% nas)  # DNA thymine
})

# ---------------------------------------------------------------------------
# get_cs_data (mocked - offline)
# ---------------------------------------------------------------------------

test_that("get_cs_data returns a data.frame with correct columns", {
  mock_df <- make_mock_stat_df()

  with_mocked_bindings(
    `.bmrb_get` = function(...) .mock_cs_stat_json,
    {
      df <- get_cs_data("ALA", "CA")
    },
    .package = "RBMRB"
  )

  expect_s3_class(df, "data.frame")
  expect_true("Val" %in% colnames(df))
  expect_true(nrow(df) >= 1L)
})

test_that("get_cs_data outlier filtering reduces row count", {
  # Build a dataset with one extreme outlier
  big_json <- .mock_cs_stat_json
  big_json$data <- c(
    big_json$data,
    list(list("99999","1","1","ALA","CA","C","9999.0","0.1","1","1","7.0","298"))
  )

  with_mocked_bindings(
    `.bmrb_get` = function(...) big_json,
    {
      df_filtered   <- get_cs_data("ALA", "CA", filtered = TRUE,  sd_limit = 2)
      df_unfiltered <- get_cs_data("ALA", "CA", filtered = FALSE)
    },
    .package = "RBMRB"
  )

  expect_true(nrow(df_filtered) < nrow(df_unfiltered))
  expect_false(9999.0 %in% df_filtered$Val)
})

test_that("get_2d_cs_data returns paired numeric vectors", {
  # Mock two separate calls for atom1 and atom2
  call_count <- 0L
  with_mocked_bindings(
    `get_cs_data` = function(residue, atom, ...) {
      call_count <<- call_count + 1L
      make_mock_stat_df()
    },
    {
      result <- get_2d_cs_data("ALA", "CA", "CB")
    },
    .package = "RBMRB"
  )

  expect_type(result$atom1, "double")
  expect_type(result$atom2, "double")
  # Lengths match (paired)
  expect_equal(length(result$atom1), length(result$atom2))
})

test_that("get_cs_from_bmrb_db with list_of_atoms delegates correctly", {
  call_log <- character(0)
  with_mocked_bindings(
    `.bmrb_get` = function(...) .mock_cs_stat_json,
    {
      df <- get_cs_from_bmrb_db(list_of_atoms = c("ALA-CB", "GLY-CA"))
    },
    .package = "RBMRB"
  )
  expect_s3_class(df, "data.frame")
})

test_that("get_cs_from_bmrb_db errors on malformed atom spec", {
  expect_error(
    get_cs_from_bmrb_db(list_of_atoms = "BADSPEC"),
    "RESIDUE-ATOM"
  )
})

# ---------------------------------------------------------------------------
# Histogram internal helpers
# ---------------------------------------------------------------------------

test_that(".resolve_atom_specs handles residue + atom", {
  specs <- RBMRB:::.resolve_atom_specs("ALA", "CA", NULL)
  expect_length(specs, 1L)
  expect_equal(specs[[1L]]$residue, "ALA")
  expect_equal(specs[[1L]]$atom,    "CA")
  expect_equal(specs[[1L]]$label,   "ALA-CA")
})

test_that(".resolve_atom_specs handles list_of_atoms", {
  specs <- RBMRB:::.resolve_atom_specs(NULL, NULL,
                                        c("ALA-CB","CYS-N","TYR-CB"))
  expect_length(specs, 3L)
  expect_equal(specs[[2L]]$residue, "CYS")
  expect_equal(specs[[2L]]$atom,    "N")
})

test_that(".resolve_atom_specs errors when nothing supplied", {
  expect_error(RBMRB:::.resolve_atom_specs(NULL, NULL, NULL))
})

test_that(".resolve_atom_specs errors on malformed list_of_atoms", {
  expect_error(RBMRB:::.resolve_atom_specs(NULL, NULL, "BADSPEC"))
})

test_that(".format_filter_desc produces readable string", {
  desc <- RBMRB:::.format_filter_desc(list(CA = 64.5, H = 7.8))
  expect_true(grepl("CA=64.50", desc))
  expect_true(grepl("H=7.80",   desc))
})

test_that(".hex_alpha converts hex colour to rgba string", {
  rgba <- RBMRB:::.hex_alpha("#1F77B4", 0.5)
  expect_match(rgba, "^rgba\\(\\d+,\\d+,\\d+,0\\.50\\)$")
})

# ---------------------------------------------------------------------------
# cs_hist (offline, mock API)
# ---------------------------------------------------------------------------

test_that("cs_hist returns invisible list with values and labels", {
  with_mocked_bindings(
    `get_cs_data` = function(...) make_mock_stat_df(),
    {
      result <- cs_hist("ALA", "CA", show_visualization = FALSE)
    },
    .package = "RBMRB"
  )
  expect_type(result, "list")
  expect_named(result, c("values","labels"))
  expect_true(length(result$values) > 0L)
  expect_equal(length(result$values), length(result$labels))
})

test_that("cs_hist with list_of_atoms works correctly", {
  with_mocked_bindings(
    `get_cs_data` = function(...) make_mock_stat_df(),
    {
      result <- cs_hist(list_of_atoms = c("ALA-CB","GLY-CA"),
                        show_visualization = FALSE)
    },
    .package = "RBMRB"
  )
  expect_true(any(result$labels == "ALA-CB"))
  expect_true(any(result$labels == "GLY-CA"))
})

test_that("cs_hist2d returns paired data", {
  with_mocked_bindings(
    `get_2d_cs_data` = function(...) list(atom1 = c(52.1, 52.5), atom2 = c(18.7, 18.9)),
    {
      result <- cs_hist2d("ALA", "CA", "CB", show_visualization = FALSE)
    },
    .package = "RBMRB"
  )
  expect_named(result, c("atom1","atom2"))
  expect_equal(length(result$atom1), 2L)
})

test_that("conditional_cs_hist returns all_values and filtered_values", {
  mock_df <- make_mock_stat_df()
  # Add Val to simulate the filtering atom data
  with_mocked_bindings(
    `get_cs_from_bmrb_db` = function(...) mock_df,
    {
      result <- conditional_cs_hist("ALA", "CB",
        filtering_rules = list(CA = 52.0),
        h_tolerance = 0.2, c_tolerance = 2.0,
        show_visualization = FALSE)
    },
    .package = "RBMRB"
  )
  expect_type(result, "list")
  expect_true("all_values"      %in% names(result))
  expect_true("filtered_values" %in% names(result))
  expect_true(length(result$all_values) >= length(result$filtered_values))
})

test_that("conditional_cs_hist errors on empty filtering_rules", {
  expect_error(
    conditional_cs_hist("ALA","CB", filtering_rules = list(),
                        show_visualization = FALSE),
    "non-empty named list"
  )
})

# ---------------------------------------------------------------------------
# Color helpers
# ---------------------------------------------------------------------------

test_that(".res_color returns a hex string for known residues", {
  col <- RBMRB:::.res_color("ALA")
  expect_match(col, "^#[0-9A-Fa-f]{6}$")
})

test_that(".res_color falls back to UNK for unknown residues", {
  col <- RBMRB:::.res_color("XYZ")
  expect_equal(unname(col), unname(RBMRB:::.AA_COLORS[["UNK"]]))
})

test_that(".axis_label returns known label for N", {
  lbl <- RBMRB:::.axis_label("N")
  expect_true(grepl("N", lbl))
})

test_that(".axis_label returns fallback for unknown atom", {
  lbl <- RBMRB:::.axis_label("WEIRD")
  expect_true(grepl("WEIRD", lbl))
})
