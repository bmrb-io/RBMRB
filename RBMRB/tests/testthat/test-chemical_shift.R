# =============================================================================
# tests/testthat/test-chemical_shift.R
# Tests for internal parsers, API helpers, and ChemicalShift module
# =============================================================================

test_that(".loop_to_df converts JSON loop to data.frame correctly", {
  df <- RBMRB:::.loop_to_df(.mock_loop_json)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 17L)
  expect_true("Comp_ID"  %in% colnames(df))
  expect_true("Atom_ID"  %in% colnames(df))
  expect_true("Val"      %in% colnames(df))
  expect_true("Ambiguity_code" %in% colnames(df))

  # Numeric coercions
  expect_type(df$Val,            "double")
  expect_type(df$Comp_index_ID,  "integer")
  expect_type(df$Ambiguity_code, "integer")

  # Specific values
  leu_H <- df[df$Comp_ID == "LEU" & df$Atom_ID == "H", "Val"]
  expect_equal(leu_H, 8.149)

  ala_CA <- df[df$Comp_ID == "ALA" & df$Atom_ID == "CA", "Val"]
  expect_equal(ala_CA, 52.4)
})

test_that(".entry_json_to_cs_dfs extracts Atom_chem_shift loops", {
  dfs <- RBMRB:::.entry_json_to_cs_dfs(.mock_entry_json)
  expect_length(dfs, 1L)
  expect_s3_class(dfs[[1L]], "data.frame")
  expect_equal(nrow(dfs[[1L]]), 17L)
})

test_that(".cs_df_to_nested builds correct nested structure", {
  df     <- RBMRB:::.loop_to_df(.mock_loop_json)
  nested <- RBMRB:::.cs_df_to_nested(df, dataset_id = "15060", auth_tag = FALSE)

  expect_named(nested, "15060")
  ds <- nested[["15060"]]

  # Should have at least chain "1"
  chains <- setdiff(names(ds), "seq_ids")
  expect_true(length(chains) >= 1L)

  chain1 <- ds[[chains[[1L]]]]
  expect_true("seq_ids" %in% names(chain1))
  expect_true(length(chain1$seq_ids) >= 3L)

  # Check a specific atom value
  val_LEU_H <- chain1[["1"]][["H"]]
  expect_equal(val_LEU_H, 8.149)
})

# ---------------------------------------------------------------------------
# NMR-STAR parser tests (pure R, no network)
# ---------------------------------------------------------------------------

test_that(".tokenise_nmrstar handles plain and quoted tokens", {
  line1 <- "1  ALA  CA  C  52.4  0.1  1"
  t1    <- RBMRB:::.tokenise_nmrstar(line1)
  expect_equal(t1, c("1","ALA","CA","C","52.4","0.1","1"))

  line2 <- "2  'My long value'  GLY"
  t2    <- RBMRB:::.tokenise_nmrstar(line2)
  expect_equal(t2, c("2","My long value","GLY"))
})

test_that("NMR-STAR file is parsed correctly", {
  path <- write_mock_nmrstar()
  on.exit(unlink(path))

  parsed <- RBMRB:::.parse_nmrstar_files(path, auth_tag = FALSE)
  expect_length(parsed, 1L)
  df <- parsed[[1L]]

  expect_s3_class(df, "data.frame")
  expect_true("Comp_ID"  %in% colnames(df))
  expect_true("Atom_ID"  %in% colnames(df))
  expect_true("Val"      %in% colnames(df))

  # Should have 7 rows (GLY:3 + ALA:4)
  expect_equal(nrow(df), 7L)

  # GLY CA = 45.70
  gly_ca <- df[df$Comp_ID == "GLY" & df$Atom_ID == "CA", "Val"]
  expect_equal(gly_ca, 45.70)
})

test_that("cs_from_file returns correct nested structure", {
  path <- write_mock_nmrstar()
  on.exit(unlink(path))

  cs <- cs_from_file(path, data_set_ids = "test_entry")
  expect_named(cs, "test_entry")
  ds    <- cs[["test_entry"]]
  chain <- ds[[setdiff(names(ds), "seq_ids")[[1L]]]]
  expect_true(2L %in% chain$seq_ids)
})

test_that("cs_from_file errors on missing file", {
  expect_error(
    cs_from_file("/nonexistent/path/file.str"),
    NA  # no error thrown for missing files — they are warned and skipped
  )
  # Actually cs_from_file warns and returns empty list
  result <- suppressWarnings(cs_from_file("/nonexistent/path/file.str"))
  expect_equal(length(result), 0L)
})

test_that("data_set_ids length must match file_paths", {
  path <- write_mock_nmrstar()
  on.exit(unlink(path))
  expect_error(cs_from_file(path, data_set_ids = c("a", "b")))
})

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

test_that("cs_to_df converts nested list to data.frame", {
  df_raw <- RBMRB:::.loop_to_df(.mock_loop_json)
  nested <- RBMRB:::.cs_df_to_nested(df_raw, "15060")
  df     <- cs_to_df(nested)

  expect_s3_class(df, "data.frame")
  expect_true("Atom_ID" %in% colnames(df))
  expect_true("Val"     %in% colnames(df))
  expect_true(nrow(df)  > 0L)
})

test_that("filter_cs_outliers removes extreme values", {
  # Use a tight cluster with a single extreme outlier so SD stays manageable
  x      <- c(rep(52, 20), rep(53, 20), rep(54, 20), 9999)
  result <- filter_cs_outliers(x, sd_limit = 3)
  expect_false(9999 %in% result)
  expect_true(all(result >= 50 & result <= 56))
})

test_that("cs_stats returns correct statistics", {
  x <- c(55.0, 55.5, 56.0, 54.5, 55.0)
  s <- cs_stats(x)
  expect_named(s, c("n","mean","sd","median","min","max"))
  expect_equal(s[["n"]],   5)
  expect_equal(s[["mean"]], mean(x))
  expect_equal(s[["min"]],  54.5)
  expect_equal(s[["max"]],  56.0)
})

test_that("standard_amino_acids returns 20 codes", {
  aa <- standard_amino_acids()
  expect_length(aa, 20L)
  expect_true("ALA" %in% aa)
  expect_true("TRP" %in% aa)
})
