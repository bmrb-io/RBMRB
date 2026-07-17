# RBMRB 3.0.0

This is a major release that modernises the package, upgrades the API
endpoint, and adds all features from PyBMRB 3.0.

## Breaking changes

* The BMRB API endpoint has moved from the deprecated
  `webapi.bmrb.wisc.edu/v2` to the current `api.bmrb.io/v2`.
  All network calls now target the new endpoint automatically.
* HTTP requests migrated from `httr` (archived) to `httr2 (>= 1.0.0)`.
* `JSON` parsing migrated from `rjson` to `jsonlite` for better null
  handling and consistent behaviour across platforms.
* Entry data is now fetched with `?format=zlib` compression, matching
  PyBMRB/pynmrstar behaviour, and decompressed with `memDecompress()`.
* The old `get_data()` function is replaced by `get_cs_data()`.
* The old `hist_cs()` function is replaced by `cs_hist()`.
* Plotting now uses `plotly` for fully interactive spectra instead of
  static `ggplot2` images.

## New features

### ChemicalShift module
* `cs_from_bmrb()` — fetch per-entry chemical shifts from BMRB as a
  nested list (dataset → chain → seq_id → atom → ppm value).
* `cs_from_file()` — parse local NMR-STAR `.str` files using a
  pure-R NMR-STAR 3 parser (no external C library required).
* `cs_table()` — flat data.frame form of per-entry chemical shifts.
* `cs_to_df()`, `cs_peek()`, `cs_summary()` — data inspection helpers.

### ChemicalShiftStatistics module
* `get_cs_data()` — database-wide chemical shift query with filtering
  by pH, temperature, ambiguity code, and outlier removal.
* `get_2d_cs_data()` — paired 2D chemical shift data for two atoms
  in the same residue.
* `get_cs_from_bmrb_db()` — raw BMRB search results (no filtering).

### Spectra module
* `n15hsqc()` — simulate an interactive 1H-15N HSQC spectrum.
* `c13hsqc()` — simulate an interactive 1H-13C HSQC spectrum.
* `tocsy()` — simulate an interactive 1H-1H TOCSY spectrum.
* `generic_2d()` — any two-atom 2D correlation with optional
  sequence-walk tracing.
* `create_n15hsqc_peaklist()`, `create_c13hsqc_peaklist()`,
  `create_tocsy_peaklist()`, `create_2d_peaklist()` — build peak
  lists without plotting.
* `export_peak_list()` — export peak lists to CSV or Sparky format.
* All spectrum functions support multi-dataset overlay, draw-trace
  between matched residues, and external peak-list overlay.

### Histogram module
* `cs_hist()` — interactive 1D chemical shift distribution as
  histogram, box plot, or violin plot.
* `cs_hist2d()` — 2D correlation heat-map or contour plot.
* `conditional_cs_hist()` — distribution conditioned on known
  chemical shifts of other atoms in the same residue.

### Utility functions
* `bmrb_test()` — three-step connectivity self-test.
* `bmrb_debug_entry()` — diagnostic for entry fetch failures.
* `bmrb_search()` — BMRB instant search.
* `bmrb_pdb_ids()` / `bmrb_from_pdb()` — cross-reference PDB IDs.
* `bmrb_citation()` — retrieve entry citations in BibTeX/text/JSON-LD.
* `bmrb_upload_entry()` — upload a local NMR-STAR file for
  temporary BMRB-hosted access (7-day token).
* `filter_cs_outliers()`, `cs_stats()` — chemical shift statistics.
* `standard_amino_acids()`, `standard_nucleic_acids()` — residue
  code reference vectors.

## Bug fixes

* Fixed silent data loss when API returns JSON `null` values in loop
  rows — `unlist()` was dropping nulls; replaced with `vapply()`.
* Fixed chemical shift search endpoint parser to use `"columns"` key
  (not `"tags"`) matching the current BMRB API response format.
* Fixed plain JSON entry response unwrapping for `{"entry_id": VALUE}`
  wrapper structure returned by the non-zlib endpoint.

---

# RBMRB 2.0.0

* Intermediate release (not on CRAN).
* Migrated API endpoint from webapi.bmrb.wisc.edu to api.bmrb.io/v2.

---

# RBMRB 1.1.0

* Added `n15hsqc()` and `c13hsqc()` functions.
* Added `get_data()` for chemical shift statistics.

---

# RBMRB 1.0.0

* Initial CRAN release.
* Basic HSQC and TOCSY simulation.
* Chemical shift histogram with `hist_cs()`.
