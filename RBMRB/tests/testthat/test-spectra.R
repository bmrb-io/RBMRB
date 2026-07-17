# =============================================================================
# tests/testthat/test-spectra.R
# Tests for peaklist builders (offline — using mock data)
# =============================================================================

# Helper: build mock dfs list as .collect_cs_dfs would return
make_mock_dfs <- function() {
  df <- make_mock_cs_df()
  list("15060" = df)
}

# ---------------------------------------------------------------------------
# N15-HSQC peaklist
# ---------------------------------------------------------------------------

test_that("create_n15hsqc_peaklist correctly pairs H and N atoms", {
  dfs <- make_mock_dfs()

  # Patch internal function to bypass network
  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_n15hsqc_peaklist(bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )

  expect_false(is.null(pl))
  expect_true(length(pl$x) >= 3L)   # GLY, ALA, LEU (or VAL, ALA, LEU depending on mock)
  expect_equal(length(pl$x), length(pl$y))
  expect_equal(length(pl$x), length(pl$dataset))
  expect_equal(length(pl$x), length(pl$residue))

  # All x values should be amide H shifts (roughly 6–11 ppm range)
  expect_true(all(pl$x > 5 & pl$x < 12))
  # All y values should be N shifts (roughly 100–135 ppm range)
  expect_true(all(pl$y > 90 & pl$y < 145))
})

test_that("create_n15hsqc_peaklist with draw_trace returns cs_track", {
  dfs1 <- make_mock_dfs()
  dfs2 <- list("15060" = make_mock_cs_df(), "99999" = make_mock_cs_df())

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs2,
    {
      pl <- create_n15hsqc_peaklist(bmrb_ids = c("15060","99999"),
                                    draw_trace = TRUE)
    },
    .package = "RBMRB"
  )
  expect_false(is.null(pl))
  expect_type(pl$cs_track, "list")
})

test_that("create_n15hsqc_peaklist returns NULL when no data", {
  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) list(),
    {
      pl <- create_n15hsqc_peaklist(bmrb_ids = "00000")
    },
    .package = "RBMRB"
  )
  expect_null(pl)
})

# ---------------------------------------------------------------------------
# C13-HSQC peaklist
# ---------------------------------------------------------------------------

test_that("create_c13hsqc_peaklist correctly pairs H and C atoms", {
  dfs <- make_mock_dfs()

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_c13hsqc_peaklist(bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )

  expect_false(is.null(pl))
  expect_equal(length(pl$x), length(pl$y))
  # H shifts: roughly 0.5–10 ppm
  expect_true(all(pl$x > 0 & pl$x < 12))
  # C shifts: roughly 10–180 ppm
  expect_true(all(pl$y > 5 & pl$y < 200))
})

test_that(".h_to_c maps hydrogen atoms to correct carbons", {
  expect_equal(RBMRB:::.h_to_c("HA"),    "CA")
  expect_equal(RBMRB:::.h_to_c("HB2"),   "CB")
  expect_equal(RBMRB:::.h_to_c("HB3"),   "CB")
  expect_equal(RBMRB:::.h_to_c("HG11"),  "CG1")
  expect_equal(RBMRB:::.h_to_c("HG12"),  "CG1")
  expect_equal(RBMRB:::.h_to_c("HG21"),  "CG2")
  expect_equal(RBMRB:::.h_to_c("HD1"),   "CD1")
  expect_equal(RBMRB:::.h_to_c("HE2"),   "CE2")
  expect_equal(RBMRB:::.h_to_c("HZ"),    "CZ")
  expect_true(is.na(RBMRB:::.h_to_c("N")))   # not a H atom
})

# ---------------------------------------------------------------------------
# TOCSY peaklist
# ---------------------------------------------------------------------------

test_that("create_tocsy_peaklist generates symmetric H-H pairs", {
  dfs <- make_mock_dfs()

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_tocsy_peaklist(bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )

  expect_false(is.null(pl))
  expect_equal(length(pl$x), length(pl$y))
  # Should include diagonal peaks (x == y) and off-diagonal (x != y)
  has_diagonal <- any(abs(pl$x - pl$y) < 1e-6)
  expect_true(has_diagonal)
  # All shifts should be proton range
  expect_true(all(pl$x > 0 & pl$x < 12))
  expect_true(all(pl$y > 0 & pl$y < 12))
})

# ---------------------------------------------------------------------------
# Generic 2D peaklist
# ---------------------------------------------------------------------------

test_that("create_2d_peaklist extracts any atom pair", {
  dfs <- make_mock_dfs()

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_2d_peaklist("N", "CA", bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )

  expect_false(is.null(pl))
  expect_equal(pl$x_atom, "N")
  expect_equal(pl$y_atom, "CA")
  expect_equal(length(pl$x), length(pl$y))

  # N shifts ~100-135, CA shifts ~40-65 for the mock data
  expect_true(all(pl$x > 90 & pl$x < 145))
  expect_true(all(pl$y > 40 & pl$y < 70))
})

test_that("create_2d_peaklist with include_next adds extra peaks", {
  dfs <- make_mock_dfs()

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl_base <- create_2d_peaklist("N", "CA", bmrb_ids = "15060",
                                     include_next = FALSE)
      pl_next <- create_2d_peaklist("N", "CA", bmrb_ids = "15060",
                                     include_next = TRUE)
    },
    .package = "RBMRB"
  )

  # With include_next we should get more peaks
  expect_true(length(pl_next$x) >= length(pl_base$x))
})

# ---------------------------------------------------------------------------
# export_peak_list
# ---------------------------------------------------------------------------

test_that("export_peak_list writes a valid CSV file", {
  dfs <- make_mock_dfs()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_n15hsqc_peaklist(bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )

  out <- export_peak_list(pl, output_file_name = tmp, output_format = "csv")
  expect_true(file.exists(tmp))
  read_back <- utils::read.csv(tmp)
  expect_equal(nrow(read_back), length(pl$x))
  expect_true("H_ppm" %in% colnames(read_back))
  expect_true("N_ppm" %in% colnames(read_back))
})

test_that("export_peak_list writes a valid Sparky file", {
  dfs <- make_mock_dfs()
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))

  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_n15hsqc_peaklist(bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )

  export_peak_list(pl, output_file_name = tmp, output_format = "sparky")
  lines <- readLines(tmp)
  expect_true(grepl("^Assignment", lines[[1L]]))
  expect_true(length(lines) > 1L)  # header + data rows
})

test_that("export_peak_list returns a data.frame invisibly", {
  dfs <- make_mock_dfs()
  with_mocked_bindings(
    `.collect_cs_dfs` = function(...) dfs,
    {
      pl <- create_n15hsqc_peaklist(bmrb_ids = "15060")
    },
    .package = "RBMRB"
  )
  out <- export_peak_list(pl)
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) > 0L)
})

# ---------------------------------------------------------------------------
# Sequence-walk builder
# ---------------------------------------------------------------------------

test_that(".build_walk produces correct segments", {
  pl <- list(
    x       = c(8.1, 7.9, 8.2, 7.8),
    y       = c(122, 119, 124, 118),
    dataset = c("A", "A", "A", "A"),
    seq_id  = c("1", "2", "3", "5"),
    residue = c("LEU","VAL","ALA","GLY")
  )
  segs <- RBMRB:::.build_walk(pl, full_walk = FALSE)
  # Gap between seq 3 and 5 → should produce 2 segments (1→2, 2→3), not 3→5
  expect_equal(length(segs), 2L)

  segs_full <- RBMRB:::.build_walk(pl, full_walk = TRUE)
  expect_equal(length(segs_full), 3L)
})
