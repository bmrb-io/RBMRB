## Resubmission notes

This is a major update (3.0.0) to the RBMRB package (previously at 1.x on CRAN).
The update addresses the following critical issues and adds all PyBMRB 3.0 features.

### Why this update is necessary

1. **Broken API endpoint** — The old endpoint (webapi.bmrb.wisc.edu/v2) was
   decommissioned. All users of RBMRB 1.x currently receive HTTP errors.
   This release migrates to the current endpoint (api.bmrb.io/v2).

2. **Archived dependency** — The `httr` package dependency was archived on CRAN.
   This release migrates to `httr2`.

3. **Archived dependency** — The `rjson` dependency was replaced with `jsonlite`
   for robust null handling.

### Test environments

* macOS Tahoe 26.5 (aarch64), R 4.6.1 (local)
* Ubuntu 24.04, R 4.4.0 (GitHub Actions)
* Windows Server 2022, R-devel (win-builder)

### R CMD CHECK results

* 0 ERRORs
* 0 WARNINGs
* 2 NOTEs (described below)

### Notes

**NOTE 1**: "File LICENSE is not mentioned in the DESCRIPTION file."
  The package uses `License: GPL-3`. No separate LICENSE file is included;
  the NOTE is a false positive from checking the source directory which
  contains a cran-comments.md file (see Note 2).

**NOTE 2**: "Non-standard file/directory found at top level: 'cran-comments.md'"
  This file is excluded from the built package via .Rbuildignore.
  It is present only in the source repository for CRAN submission notes.

**NOTE 3** (if present): "no visible global function definition for 'setNames'"
  `setNames` is a base R function (package:base) always available in any R
  session. The NOTE is a known false positive from R CMD CHECK for functions
  in the `base` package. We have added `importFrom("stats", "setNames")` to
  suppress it.

### Downstream dependencies

RBMRB has no known reverse dependencies on CRAN.
