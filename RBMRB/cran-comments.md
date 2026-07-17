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


### R CMD CHECK results

* 0 ERRORs
* 0 WARNINGs
* 0 NOTEs 

### Downstream dependencies

RBMRB has no known reverse dependencies on CRAN.
