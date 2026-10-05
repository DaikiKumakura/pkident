# pkident

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/pkident/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/pkident/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Structural and practical identifiability for pharmacokinetic and pharmacodynamic ODE models, built on rxode2.

pkident answers three questions about a model before (or after) it is fitted:

1. **Structure.** With continuous, noise-free observation of the chosen outputs after the known doses, which parameters can be determined at all? Which combinations (for example `V/F`) can?
2. **Design.** With these sampling times, doses and assay error, how precisely would each parameter be estimated? Which added measurement, sampling time or dose level would help?
3. **Fit.** In an nlmixr2 fit, which parameters do the data actually determine, and is a flat profile a property of the model or of the data?

Every decision is **local** and says so, and anything the methods cannot decide is reported as `unresolved` rather than guessed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("DaikiKumakura/pkident")
```

pkident needs rxode2 and symengine. `classify_profiles()` also needs nlmixr2est and nlmixr2extra.

## Structural identifiability

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

Two independent methods must agree: Lie derivatives of the outputs evaluated at the known initial state (which keeps the information a dose gives about volumes), and the rank and null space of output sensitivities solved with rxode2. Known input signals, such as a plasma concentration driving a pharmacodynamic model, are supported.

## Practical identifiability of a design

```r
practical_identifiability(
  reference_model("A3"),                       # oral, F known
  design = c(0.5, 1, 2),
  values = c(ka = 1, CL = 2, V = 20, DOSE = 100),
  error  = residual_error(add = 0.05, prop = 0.1)
)
#> Condition number (log scale): 36900
#>   CL           RSE    975.1%  unresolved
#>   V            RSE    224.2%  unresolved
#>   ka           RSE    257.2%  unresolved
#> Weakest direction: +0.94 log(CL) -0.25 log(ka) -0.22 log(V) (expected SD on log scale 10.3)
```

Three samples in the first two hours do not observe elimination; eight samples up to 24 hours give RSEs of 5–13%. `compare_designs()` puts candidate designs side by side, with the structural status for the outputs each candidate observes. Doses can be given as event tables (repeated infusions), and residual error can be additive, proportional, combined or exponential.

## Likelihood profiles of an nlmixr2 fit

```r
classify_profiles(fit, structural = structural_identifiability(emax_model),
                  map = c(te0 = "E0", temax = "Emax", tec50 = "EC50"))   # excerpt
#>   temax  non_identifiable  open towards upper values (explored -0.298 to 3.7, 10 fits)
#>          -> structurally identifiable but not determined by these data (structural: identifiable)
```

Profiles are classified as bounded, one-sided, flat or shallow (Raue et al. 2009). Each profile point is re-estimated from a neighbouring point when needed, because a stuck re-estimation can only overstate the objective function.

## Validation

| Check | Result |
| --- | --- |
| Ten reference models with expected results fixed from the literature or by derivation (and one variant) | All decisions as expected at four seeds; the two structural methods agree at every point |
| Fisher information against the analytic result | Relative difference 1e-11 (also for repeated doses from an event table) |
| Expected RSEs against 1,000 repeated estimates on simulated data | Ratio 0.96–1.08 in four scenarios; limitations at RSEs of 30–40% documented |
| Profile classification on three synthetic nlmixr2 fits | All ten parameters classified as expected |

Details: `validation/` in this repository.

## Vignettes and case studies

* *Structural identifiability of reference PK/PD models*
* *Case study: nimotuzumab target-mediated disposition* — with the public `nimoData` (nlmixr2data), clearance, volume and koff are determined, kon, kdeg, target baseline and kint are not, and measuring total target would determine six of the seven parameters.
* *Case study: PSA decline and regrowth* — the regrowth rate becomes estimable only about twice the time to nadir after the start of treatment.

The analysis scripts are in `case-studies/`.

## Scope and limits

* Decisions are local. A locally identifiable parameter can have a finite number of alternative values (for example the flip-flop of absorption and elimination).
* Fixed effects of a typical individual. Identifiability of between-subject variability and population designs are not covered (see PopED for optimal design).
* Closed-form models (`linCmt()`) must be written as ODEs.
* The expected information is a local, linear approximation; it is reliable when RSEs are below about 30%, and profiles should be examined otherwise.

## License

GPL-3
