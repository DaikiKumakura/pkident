# pkident

<!-- badges: start -->
[![R-universe](https://daikikumakura.r-universe.dev/badges/pkident)](https://daikikumakura.r-universe.dev/pkident)
[![R-CMD-check](https://github.com/DaikiKumakura/pkident/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/pkident/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210471.svg)](https://doi.org/10.5281/zenodo.24210471)
<!-- badges: end -->

**Which parameters of a PK/PD model can the data determine? Structural and
practical identifiability for ODE models written in rxode2.**

pkident checks a pharmacokinetic or pharmacodynamic model before or after it
is fitted:

- **Structure:** with the chosen outputs and known doses, which parameters
  can be determined at all, and which combinations (such as `CL/F`) can?
- **Design:** with these sampling times, doses and assay error, how precisely
  would each parameter be estimated, and which change to the design helps?
- **Fit:** in an nlmixr2 fit, which parameters do the data actually determine,
  and is a flat profile a property of the model or of the data?

Every decision is **local** and says so. Anything the methods cannot decide
is reported as `unresolved` rather than guessed, and two independent methods
must agree before a structural decision is made.

**Documentation:** <https://daikikumakura.github.io/pkident/> (function
reference and three case studies)

## Workflow

```text
   model (rxode2 ODEs) + outputs + doses
                 │
   pkpd_model()                       describe what is observed and known
                 │
   structural_identifiability()       can the parameters be determined at all?
                 │
   practical_identifiability()        how precisely, with this design?
   compare_designs()                  which design is better?
                 │
   classify_profiles()                what did an nlmixr2 fit determine?
```

| Question | Function | Answer | Needs |
|---|---|---|---|
| Describe the model | `pkpd_model()`, `bolus_dose()`, `infusion_dose()`, `reference_model()` | a model object with outputs, doses and known quantities | rxode2 |
| Structure | `structural_identifiability()` | `identifiable`, `combination_only`, `non_identifiable` or `unresolved` per parameter, with the identifiable combinations | symengine |
| Design | `practical_identifiability()`, `compare_designs()`, `experiment()`, `residual_error()` | expected RSE per parameter, condition number, weakest direction | — |
| Fit | `classify_profiles()` | bounded, one-sided, flat or shallow profile per parameter, linked to the structural result | nlmixr2est, nlmixr2extra |
| Report | `as_table()` | a data frame of any result | — |

## Installation

```r
install.packages("pkident", repos = c(
  "https://daikikumakura.r-universe.dev",
  "https://cloud.r-project.org"
))
```

Or from GitHub: `remotes::install_github("DaikiKumakura/pkident")`. pkident
needs rxode2 and symengine; `classify_profiles()` also needs nlmixr2est and
nlmixr2extra.

## Quick start: structure

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

Only `ka` and the ratios `CL/F` and `CL/V` are determined. The two methods
are Lie derivatives of the outputs at the known initial state (which keeps
what a dose tells about volumes) and the rank and null space of output
sensitivities solved with rxode2. Known input signals, such as a plasma
concentration driving a PD model, are supported.

## Design: how precise would the estimates be?

The same model with `F` known (`reference_model("A3")`), sampled three times
in the first two hours, and then eight times up to 24 hours:

```r
values <- c(ka = 1, CL = 2, V = 20, DOSE = 100)
error  <- residual_error(add = 0.05, prop = 0.1)

practical_identifiability(reference_model("A3"), design = c(0.5, 1, 2),
                          values = values, error = error)
#> Condition number (log scale): 36900
#>   CL           RSE    975.1%  unresolved
#>   V            RSE    224.2%  unresolved
#>   ka           RSE    257.2%  unresolved
#> Weakest direction: +0.94 log(CL) -0.25 log(ka) -0.22 log(V) (expected SD on log scale 10.3)

practical_identifiability(reference_model("A3"), design = c(0.5, 1, 2, 4, 6, 8, 12, 24),
                          values = values, error = error)
#> Condition number (log scale): 13.3
#>   CL           RSE      4.7%  identifiable
#>   V            RSE      7.1%  identifiable
#>   ka           RSE     14.6%  identifiable
```

Three early samples do not observe elimination. `compare_designs()` puts
candidate designs side by side, with the structural status for the outputs
each one observes. Doses can be event tables (for example repeated
infusions), and residual error can be additive, proportional, combined or
exponential.

## Fit: what did the data determine?

```r
classify_profiles(fit, structural = structural_identifiability(emax_model),
                  map = c(te0 = "E0", temax = "Emax", tec50 = "EC50"))   # excerpt
#>   temax  non_identifiable  open towards upper values (explored -0.298 to 3.7, 10 fits)
#>          -> structurally identifiable but not determined by these data (structural: identifiable)
```

Profiles are classified as bounded, one-sided, flat or shallow (Raue et al.
2009) and compared with the structural result, which tells a model problem
from a data problem. A profile point is re-estimated from a neighbouring point
when needed, because a stuck re-estimation can only overstate the objective
function.

## Validation

| Check | Result |
| --- | --- |
| Ten reference models with expected results fixed from the literature or by derivation (and one variant) | All decisions as expected at four seeds; the two structural methods agree at every point |
| Fisher information against the analytic result | Relative difference 1e-11 (also for repeated doses from an event table) |
| Expected RSEs against 1,000 repeated estimates on simulated data | Ratio 0.96–1.08 in four scenarios; limitations at RSEs of 30–40% documented |
| Profile classification on three synthetic nlmixr2 fits | All ten parameters classified as expected |

Details are in [`validation/`](https://github.com/DaikiKumakura/pkident/tree/master/validation). R CMD check passes in CI.

## Case studies

- *Structural identifiability of reference PK/PD models.*
- *Nimotuzumab target-mediated disposition:* with the public `nimoData`
  (nlmixr2data), clearance, volume and koff are determined, kon, kdeg, target
  baseline and kint are not, and measuring total target would determine six
  of the seven parameters.
- *PSA decline and regrowth:* the regrowth rate becomes estimable only about
  twice the time to nadir after the start of treatment.

They are articles on the documentation site; the analysis scripts are in
[`case-studies/`](https://github.com/DaikiKumakura/pkident/tree/master/case-studies).

## Scope and limits

- Decisions are local. A locally identifiable parameter can have a finite
  number of alternative values (for example the flip-flop of absorption and
  elimination).
- Fixed effects of a typical individual. Identifiability of between-subject
  variability and population designs are not covered (see PopED for optimal
  design).
- Closed-form models (`linCmt()`) must be written as ODEs.
- The expected information is a local, linear approximation: reliable when
  RSEs are below about 30%; examine likelihood profiles otherwise.

## Related reading

- Glossary: [structural identifiability](https://daikikumakura.github.io/glossary/structural-identifiability.html),
  [practical identifiability](https://daikikumakura.github.io/glossary/practical-identifiability.html),
  [profile likelihood](https://daikikumakura.github.io/glossary/profile-likelihood.html).
- Analyses on what PK/PD data can identify:
  [warfarin PK/PD reanalysis](https://daikikumakura.github.io/writing/warfarin-pkpd.html),
  [nimotuzumab population PK](https://daikikumakura.github.io/writing/nimotuzumab-poppk-en.html),
  [tarlatamab summary PK](https://daikikumakura.github.io/writing/tarlatamab-summary-pk-en.html).

## Citation

Cite the archived version: <https://doi.org/10.5281/zenodo.24210471>. The
citation metadata is in [`CITATION.cff`](https://github.com/DaikiKumakura/pkident/blob/master/CITATION.cff).

## License

GPL-3
