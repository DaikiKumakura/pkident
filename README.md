# pkident

<!-- badges: start -->
[![R-universe](https://daikikumakura.r-universe.dev/badges/pkident)](https://daikikumakura.r-universe.dev/pkident)
[![R-CMD-check](https://github.com/DaikiKumakura/pkident/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/pkident/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210471.svg)](https://doi.org/10.5281/zenodo.24210471)
<!-- badges: end -->

**Which parameters of a PK/PD model can the data determine? Structural and
practical identifiability for ODE models written in rxode2.**

## Installation

```r
install.packages("pkident", repos = c(
  "https://daikikumakura.r-universe.dev",
  "https://cloud.r-project.org"
))
```

<a id="quick-start-structure"></a>

## Quick start

An oral one-compartment model with unknown bioavailability:

```r
library(pkident)

m <- pkpd_model(
  "d/dt(depot)   = -ka*depot
   d/dt(central) = ka*depot - CL/V*central
   cp            = central/V",
  outputs = "cp",
  doses   = bolus_dose("depot", "F*DOSE"),
  known   = "DOSE"
)

structural_identifiability(m)
#> <pkident structural identifiability>
#> Decision is LOCAL (Lie derivatives at the known initial state up to order 8, checked by output sensitivities; 5 random points).
#> Methods agree at 5 of 5 points.
#>   CL           combination_only
#>   F            combination_only
#>   V            combination_only
#>   ka           identifiable
#> Identifiable combinations:
#>   CL/F
#>   CL/V
```

## Documentation

[Package website](https://daikikumakura.github.io/pkident/) · [Function reference](https://daikikumakura.github.io/pkident/reference/) · [Changelog](https://daikikumakura.github.io/pkident/news/) · [Releases](https://github.com/DaikiKumakura/pkident/releases) · [R-universe](https://daikikumakura.r-universe.dev/pkident) · [Usage guide](https://github.com/DaikiKumakura/pkident/blob/master/README-GUIDE.md)

## Validation and limitations

- **Validation:** Checked identifiability classifications against reference models and compared Fisher information with analytical and simulated results.
- **Results:** Both structural methods matched the expected classifications; the analytical Fisher information comparison agreed to 1e-11 relative error.
- **Limitations:** Design results depend on the model, parameters and sampling assumptions; local uncertainty does not establish global identifiability.
- **Reproduction:** [vignettes/reference-models.Rmd](https://github.com/DaikiKumakura/pkident/blob/master/vignettes/reference-models.Rmd).

## Citation

Kumakura, D. (2026). *pkident: Structural and Practical Identifiability for Pharmacokinetic and Pharmacodynamic Models* [Software]. [DOI: 10.5281/zenodo.24210471](https://doi.org/10.5281/zenodo.24210471).

Citation metadata: [CITATION.cff](CITATION.cff).

## License

This project is licensed under GPL-3.0-only. See [LICENSE.md](LICENSE.md) for the license terms.
