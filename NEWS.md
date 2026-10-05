# pkident 0.1.0

First release. Decisions are local and are stated as such in every result.

## Model definition

* `pkpd_model()` reads ODE models written for rxode2 (code, compiled models or
  model functions) and records outputs, doses (`bolus_dose()`,
  `infusion_dose()`), known input signals, known constants and initial
  conditions.
* `reference_models()` and `reference_model()` provide ten acceptance models
  with expected results fixed from primary literature or analytical
  derivation.

## Structural identifiability

* `structural_identifiability()` decides local structural identifiability with
  Lie derivatives evaluated at the known initial state, checked by the rank and
  null space of output sensitivities solved with rxode2. A decision is made
  only if all random points and both methods agree; otherwise the result is
  `unresolved`. Identifiable monomial combinations (for example `V/F`) are
  reported.

## Practical identifiability and design

* `practical_identifiability()` computes the expected Fisher information for a
  sampling design with additive, proportional, combined or exponential residual
  error (`residual_error()`), with doses from the model or from an event table.
* `compare_designs()` and `experiment()` compare user-supplied candidate
  designs: added outputs, sampling times or dose levels.
* `classify_profiles()` classifies likelihood profiles of nlmixr2 fits
  (bounded, one-sided, flat, shallow; Raue et al. 2009) using
  `nlmixr2extra::profileFixed()` with warm-started re-estimation, and sets
  them against the structural result.

## Validation

* All ten acceptance models (and a variant) agree with their expected results
  at four seeds, and the two structural methods agree at every point.
* Expected RSEs agree with the spread of 1,000 repeated estimates within
  0.96–1.08 in four scenarios; limitations at RSEs of 30–40% are documented.
* Profile classification is correct for all ten parameters of three synthetic
  nlmixr2 fits.
* Vignettes: reference models, nimotuzumab TMDD (public `nimoData`), PSA
  decline and regrowth.
