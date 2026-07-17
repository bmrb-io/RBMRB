# =============================================================================
# tests/testthat/test-api_utils.R
# Tests for public API utility functions (all mocked for offline running)
# =============================================================================

test_that("bmrb_list_entries returns a character vector", {
  with_mocked_bindings(
    `.bmrb_get` = function(...) list("10001","10002","10003"),
    {
      ids <- bmrb_list_entries("macromolecules")
    },
    .package = "RBMRB"
  )
  expect_type(ids, "character")
  expect_true(length(ids) >= 1L)
})

test_that("bmrb_list_entries validates database argument", {
  expect_error(bmrb_list_entries("invalid_db"))
})

test_that("bmrb_search returns a data.frame", {
  mock_search <- list(
    list(value = "15060", label = "Ubiquitin", sub_date = "2006-09-07",
         authors = list("Smith", "Jones"), link = "/data/entry/15060"),
    list(value = "15061", label = "Lysozyme",  sub_date = "2007-01-01",
         authors = list("Brown"),             link = "/data/entry/15061")
  )
  with_mocked_bindings(
    `.bmrb_get` = function(...) mock_search,
    {
      df <- bmrb_search("ubiquitin")
    },
    .package = "RBMRB"
  )
  expect_s3_class(df, "data.frame")
  expect_true("value"    %in% colnames(df))
  expect_true("label"    %in% colnames(df))
  expect_true("sub_date" %in% colnames(df))
  expect_equal(nrow(df), 2L)
})

test_that("bmrb_search returns empty data.frame for no results", {
  with_mocked_bindings(
    `.bmrb_get` = function(...) list(),
    {
      df <- bmrb_search("zzznotfound")
    },
    .package = "RBMRB"
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
})

test_that("bmrb_pdb_ids returns a data.frame with correct columns", {
  mock_pdb <- list(
    list(pdb_id = "2JM0", match_type = "Exact",      comment = NULL),
    list(pdb_id = "1XYZ", match_type = "BLAST Match", comment = "similar")
  )
  with_mocked_bindings(
    `.bmrb_get` = function(...) mock_pdb,
    {
      df <- bmrb_pdb_ids(15000)
    },
    .package = "RBMRB"
  )
  expect_s3_class(df, "data.frame")
  expect_named(df, c("pdb_id","match_type","comment"))
  expect_equal(nrow(df), 2L)
  expect_equal(df$pdb_id[[1L]], "2JM0")
})

test_that("bmrb_pdb_ids returns empty data.frame when no PDB links", {
  with_mocked_bindings(
    `.bmrb_get` = function(...) list(),
    {
      df <- bmrb_pdb_ids(99999)
    },
    .package = "RBMRB"
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
})

test_that("bmrb_from_pdb returns a data.frame with correct columns", {
  mock_bmrb <- list(
    list(bmrb_id = "15000", match_type = "Exact", comment = NULL)
  )
  with_mocked_bindings(
    `.bmrb_get` = function(...) mock_bmrb,
    {
      df <- bmrb_from_pdb("2JM0")
    },
    .package = "RBMRB"
  )
  expect_s3_class(df, "data.frame")
  expect_named(df, c("bmrb_id","match_type","comment"))
  expect_equal(df$bmrb_id[[1L]], "15000")
})

test_that(".null_na handles NULL and non-NULL correctly", {
  expect_equal(RBMRB:::.null_na(NULL),  NA_character_)
  expect_equal(RBMRB:::.null_na("abc"), "abc")
  expect_equal(RBMRB:::.null_na(123),   "123")
})

test_that(".ensure_vector flattens single-element list", {
  expect_equal(RBMRB:::.ensure_vector(list("a")), "a")
  expect_equal(RBMRB:::.ensure_vector(c(1L, 2L, 3L)), c(1L, 2L, 3L))
})

test_that(".is_metabolomics detects bmse prefix", {
  expect_true(RBMRB:::.is_metabolomics("bmse000034"))
  expect_true(RBMRB:::.is_metabolomics("BMSE000034"))
  expect_false(RBMRB:::.is_metabolomics("15060"))
  expect_false(RBMRB:::.is_metabolomics(15060))
})

test_that("bmrb_status returns a list", {
  mock_status <- list(
    server_version = "2.0", databases = list(
      macromolecules = list(num_entries = 12000L, num_chemical_shifts = 5e6)
    )
  )
  with_mocked_bindings(
    `.bmrb_get` = function(...) mock_status,
    {
      s <- bmrb_status()
    },
    .package = "RBMRB"
  )
  expect_type(s, "list")
  expect_true("databases" %in% names(s))
})

test_that("bmrb_upload_entry errors when file does not exist", {
  expect_error(bmrb_upload_entry("/nonexistent/file.str"), "File not found")
})

test_that("bmrb_fetch_entry returns a named list", {
  with_mocked_bindings(
    `.bmrb_get` = function(...) .mock_entry_json,
    {
      result <- bmrb_fetch_entry("15060")
    },
    .package = "RBMRB"
  )
  expect_type(result, "list")
  expect_named(result, "15060")
})
