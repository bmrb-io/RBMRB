# =============================================================================
# tests/testthat/helper-mocks.R
# Shared mock data and HTTP interceptors for rBMRB tests
# =============================================================================

# ---------------------------------------------------------------------------
# Minimal fake Atom_chem_shift data (entry 15060-like)
# ---------------------------------------------------------------------------

.mock_atom_cs_rows <- list(
  list("1","1","1","1","1","1","LEU","H", "H","1","8.149","0.005",".", "1",".",".",".", ".","20","LEU","HN",".", "15060","1"),
  list("2","1","1","1","1","1","LEU","CA","C","13","56.016","0.10",".", "1",".",".",".", ".","20","LEU","CA",".", "15060","1"),
  list("3","1","1","1","1","1","LEU","CB","C","13","42.180","0.10",".", "1",".",".",".", ".","20","LEU","CB",".", "15060","1"),
  list("4","1","1","1","1","1","LEU","N", "N","15","122.739","0.10",".","1",".",".",".", ".","20","LEU","N", ".", "15060","1"),
  list("5","1","1","2","2","2","VAL","H", "H","1","8.048","0.005",".", "1",".",".",".", ".","21","VAL","HN",".", "15060","1"),
  list("6","1","1","2","2","2","VAL","CA","C","13","63.412","0.10",".", "1",".",".",".", ".","21","VAL","CA",".", "15060","1"),
  list("7","1","1","2","2","2","VAL","CB","C","13","32.100","0.10",".", "1",".",".",".", ".","21","VAL","CB",".", "15060","1"),
  list("8","1","1","2","2","2","VAL","N", "N","15","119.200","0.10",".","1",".",".",".", ".","21","VAL","N", ".", "15060","1"),
  list("9","1","1","3","3","3","ALA","H", "H","1","7.900","0.005",".", "1",".",".",".", ".","22","ALA","HN",".", "15060","1"),
  list("10","1","1","3","3","3","ALA","CA","C","13","52.400","0.10",".", "1",".",".",".", ".","22","ALA","CA",".", "15060","1"),
  list("11","1","1","3","3","3","ALA","CB","C","13","18.700","0.10",".", "1",".",".",".", ".","22","ALA","CB",".", "15060","1"),
  list("12","1","1","3","3","3","ALA","N", "N","15","124.100","0.10",".","1",".",".",".", ".","22","ALA","N", ".", "15060","1"),
  list("13","1","1","1","1","1","LEU","HA","H","1","4.215","0.005",".","1",".",".",".", ".","20","LEU","HA",".", "15060","1"),
  list("14","1","1","1","1","1","LEU","HB2","H","1","1.700","0.01",".","1",".",".",".", ".","20","LEU","HB2",".", "15060","1"),
  list("15","1","1","2","2","2","VAL","HA","H","1","3.980","0.005",".","1",".",".",".", ".","21","VAL","HA",".", "15060","1"),
  list("16","1","1","3","3","3","ALA","HA","H","1","4.120","0.005",".","1",".",".",".", ".","22","ALA","HA",".", "15060","1"),
  list("17","1","1","3","3","3","ALA","HB","H","1","1.390","0.01",".","1",".",".",".", ".","22","ALA","HB",".", "15060","1")
)

.mock_cs_tags <- c(
  "ID","Assembly_atom_ID","Entity_assembly_ID","Entity_ID",
  "Comp_index_ID","Seq_ID","Comp_ID","Atom_ID","Atom_type",
  "Atom_isotope_number","Val","Val_err","Assign_fig_of_merit",
  "Ambiguity_code","Occupancy","Resonance_ID",
  "Auth_entity_assembly_ID","Auth_asym_ID","Auth_seq_ID",
  "Auth_comp_ID","Auth_atom_ID","Details","Entry_ID",
  "Assigned_chem_shift_list_ID"
)

.mock_loop_json <- list(
  category = "_Atom_chem_shift",
  tags     = as.list(.mock_cs_tags),
  data     = .mock_atom_cs_rows
)

.mock_entry_json <- list(
  bmrb_id    = "15060",
  saveframes = list(
    list(
      category = "assigned_chemical_shifts",
      loops    = list(.mock_loop_json)
    )
  )
)

# ---------------------------------------------------------------------------
# Minimal NMR-STAR file text for parser tests
# ---------------------------------------------------------------------------

.mock_nmrstar_text <- c(
  "data_test",
  "save_frame_1",
  "_Entry.ID 99999",
  "loop_",
  "_Atom_chem_shift.ID",
  "_Atom_chem_shift.Comp_index_ID",
  "_Atom_chem_shift.Comp_ID",
  "_Atom_chem_shift.Atom_ID",
  "_Atom_chem_shift.Atom_type",
  "_Atom_chem_shift.Val",
  "_Atom_chem_shift.Val_err",
  "_Atom_chem_shift.Ambiguity_code",
  "_Atom_chem_shift.Entity_assembly_ID",
  "1  1  GLY  H   H  8.20  0.01  1  1",
  "2  1  GLY  CA  C  45.70 0.10  1  1",
  "3  1  GLY  N   N  109.3 0.10  1  1",
  "4  2  ALA  H   H  7.85  0.01  1  1",
  "5  2  ALA  CA  C  52.40 0.10  1  1",
  "6  2  ALA  CB  C  18.70 0.10  1  1",
  "7  2  ALA  N   N  124.1 0.10  1  1",
  "stop_",
  "save_"
)

# ---------------------------------------------------------------------------
# Minimal CS statistics mock response
# ---------------------------------------------------------------------------

.mock_cs_stat_tags <- c("Entry_ID","Entity_assembly_ID","Comp_index_ID",
                         "Comp_ID","Atom_ID","Atom_type","Val","Val_err",
                         "Ambiguity_code","Assigned_chem_shift_list_ID",
                         "pH","Temperature")

.mock_cs_stat_rows <- list(
  list("10001","1","1","ALA","CA","C","52.1","0.1","1","1","7.0","298"),
  list("10002","1","1","ALA","CA","C","52.5","0.1","1","1","7.0","298"),
  list("10003","1","1","ALA","CA","C","51.9","0.1","1","1","6.8","310"),
  list("10004","1","1","ALA","CA","C","53.0","0.1","1","1","7.2","298"),
  list("10005","1","1","ALA","CA","C","52.3","0.1","1","1","7.0","298")
)

.mock_cs_stat_json <- list(
  tags = as.list(.mock_cs_stat_tags),
  data = .mock_cs_stat_rows
)

# ---------------------------------------------------------------------------
# Helper to build a mock data.frame directly (bypasses API)
# ---------------------------------------------------------------------------

make_mock_cs_df <- function() {
  RBMRB:::.loop_to_df(.mock_loop_json)
}

make_mock_stat_df <- function() {
  RBMRB:::.loop_to_df(.mock_cs_stat_json)
}

# ---------------------------------------------------------------------------
# Write a temporary NMR-STAR file
# ---------------------------------------------------------------------------

write_mock_nmrstar <- function(dir = tempdir()) {
  path <- file.path(dir, "mock_entry.str")
  writeLines(.mock_nmrstar_text, path)
  path
}
