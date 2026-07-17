# RBMRB <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/RBMRB)](https://CRAN.R-project.org/package=RBMRB)
<!-- badges: end -->

**RBMRB** is a comprehensive R interface to the
[Biological Magnetic Resonance Data Bank (BMRB)](https://bmrb.io),
mirroring all functionality of [PyBMRB 3.0](https://github.com/bmrb-io/PyBMRB).

## Installation

```r
# From CRAN
install.packages("RBMRB")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("uwbmrb/RBMRB")
```

## Quick start

```r
library(RBMRB)

# Verify connectivity
bmrb_test()

# Fetch chemical shifts
cs <- cs_from_bmrb(15060)
cs_peek(cs)

# Simulate spectra (interactive Plotly)
n15hsqc(bmrb_ids = 15060, legend = "residue")
c13hsqc(bmrb_ids = 15060)
tocsy(bmrb_ids   = 15060)

# Compare two entries
n15hsqc(bmrb_ids = c(17076, 17077), draw_trace = TRUE)

# Chemical shift histograms
cs_hist("ALA", "CA")
cs_hist2d("ALA", "CA", "CB")
conditional_cs_hist("SER", "N", filtering_rules = list(CA = 55.5))

# Export peak list (CSV or Sparky)
pl <- create_n15hsqc_peaklist(bmrb_ids = 15060)
export_peak_list(pl, "peaks.csv")
export_peak_list(pl, "peaks.sparky", output_format = "sparky")
```

## Modules

| Module | Functions | Description |
|--------|-----------|-------------|
| **ChemicalShift** | `cs_from_bmrb()`, `cs_from_file()` | Per-entry chemical shifts |
| **ChemicalShiftStatistics** | `get_cs_data()`, `get_2d_cs_data()` | Database-wide statistics |
| **Spectra** | `n15hsqc()`, `c13hsqc()`, `tocsy()`, `generic_2d()` | Spectrum simulation |
| **Histogram** | `cs_hist()`, `cs_hist2d()`, `conditional_cs_hist()` | Distributions |

## Citation

If you use RBMRB in published research, please cite:

> Ulrich EL et al. (2008). BioMagResBank. *Nucleic Acids Research*
> **36**, D402–D408. https://doi.org/10.1093/nar/gkm957

> Baskaran K et al. (2021). PyBMRB: Data visualization tool for
> BioMagResBank. *Proc. 20th Python in Science Conf.*, 59–62.
