# pkident — Specification v0.1

Status: specification agreed on 2026-10-04. Repository not yet created. It becomes public once milestone S3 (structural identifiability) works.

## 1. Purpose

An R package that determines, for a pharmacokinetic/pharmacodynamic ODE model, **which parameters (or combinations of parameters) can be determined from a given dosing and observation design**, and, when some cannot, **which additional measurement would make them determinable**.

The decision it supports: *before estimating anything, can the parameters we care about be identified with this model structure and this study design? If not, should we simplify the structure or add a measurement?* The package turns the habit "check identifiability before adding parameters" into a tool.

Intended use:

```r
library(pkident)

m <- pkpd_model(model, outputs = "cp", doses = bolus_dose("depot"))

structural_identifiability(m)            # parameters not determinable even with perfect data
practical_identifiability(m, design)     # how weakly they are determined with this sampling and error
compare_designs(m, candidates)           # which added measurement would make them determinable
classify_profiles(fit)                   # classify likelihood profiles of an actual fit (uses nlmixr2extra)
```

## 2. Scope

**In scope (v0.1)**

- ODE models written in rxode2/nlmixr2 syntax (`d/dt()`), with observation equations that are functions of states and parameters (for example `cp = central / v`). Several outputs are allowed.
- Known doses: instantaneous bolus (initial condition), zero-order infusion (known input), and first-order absorption from a depot.
- Known input signals, such as a plasma concentration `Cp(t)` driving a pharmacodynamic model. Their time derivatives are treated as known (added after milestone S1; required by the pharmacodynamic reference models).
- Fixed effects of a typical individual.

**Out of scope (v0.1)**

- Estimation, model selection, covariate selection, clinical interpretation and dose recommendation.
- Optimal design (the domain of PopED). The package only compares candidates supplied by the user.
- Computing likelihood profiles (done by `nlmixr2extra::profile()`).
- Structural identifiability of between-subject variability. Considered for v0.2.
- Closed-form models (`linCmt()`). v0.1 accepts models rewritten as ODEs.
- Proof of global identifiability. v0.1 decides **local identifiability only** and states this in every result.

## 3. Position relative to existing tools (checked 2026-10-04)

| Tool | What it does | Relation |
| --- | --- | --- |
| nlmixr2extra (`profile()`, `profileLlp()`, `profileFixed()`) | Likelihood profiling and profile confidence intervals for nlmixr2 FOCEi fits | Used for computation. This package classifies the profiles and reconciles them with the structural result |
| FME (CRAN) | Sensitivity and collinearity for generic ODE models | Not designed around doses, PK/PD outputs or mixed-effects models |
| StructuralIdentifiability.jl, SIAN (Julia); STRIKE-GOLDD (MATLAB) | Structural identifiability | Not available in R or in a PK/PD workflow. Not used for validation |
| PopED (CRAN) | Population Fisher information and optimal design | No optimization here; diagnosis of identifiability only |
| Pharmpy | Pharmacometrics infrastructure (Python) | No identifiability functionality (code search) |

## 4. Definitions

| Term | Definition |
| --- | --- |
| Local structural identifiability | With noise-free continuous observation, parameters are uniquely determined in a neighbourhood of the true values |
| Global identifiability | Uniquely determined over the whole parameter space. Not decided in v0.1 (for example, flip-flop in the one-compartment oral model is locally identifiable but has two global solutions) |
| Practical identifiability | Parameters are determined with finite confidence intervals given finite sampling and measurement error |
| Identifiable combination | A function of parameters (for example `V/F`) that is determined although its members are not |

Every result uses only these statuses and never fills gaps by guessing:

| Status | Meaning |
| --- | --- |
| `identifiable` | Locally identifiable under the stated method |
| `non_identifiable` | Not identifiable; the parameters that move together are reported |
| `combination_only` | Not identifiable alone, but identifiable within the reported combination |
| `unresolved` | Could not be decided (computation limit reached, poor numerical conditioning, or the two methods disagree) |

## 5. Methods

### 5.1 Structural identifiability (core)

**Primary method: Lie derivatives of the outputs evaluated at the known initial state (Taylor-series approach).** Revised at milestone S3 (see below).

1. Compute Lie derivatives of the outputs symbolically with `symengine`, using the symbolic representation available from `rxode2`. Known input signals contribute their own time derivatives as additional known symbols.
2. Substitute the known initial state (set by doses and initial conditions) at time zero. The output derivatives at time zero are then functions of the parameters only.
3. Evaluate the Jacobian of these output derivatives with respect to the parameters, on a log scale and with normalized rows, at several random points, and determine its rank. Increase the derivative order until the rank is full, or until it has not increased for three consecutive orders (rank deficiency), or until a limit is reached (`unresolved`).
4. If the rank is deficient, use the null space of the Jacobian to report the directions along which parameters can move without changing the output, the parameters that are not identifiable alone, and identifiable monomial combinations (for example `V/F`).
5. Known input signals are given random values for themselves and their derivatives (a generic input).

Why the initial state is used: in pharmacokinetic models the dose fixes the initial state, and this is often what makes a volume identifiable (for example `y(0) = DOSE/V`). The generic-state rank used by STRIKE-GOLDD treats initial states as unknown and would miss this information. Rank deficiency at a finite derivative order is a strong but not a formal proof of non-identifiability; the plateau rule and the agreement across points are recorded with every result, and the secondary method (S4) is a further check.

**Secondary method: numerical check from sensitivities.** Solve the sensitivity equations with rxode2 and determine the rank of the noise-free output sensitivity matrix on a dense time grid at the same random points. If the two methods disagree, the result is `unresolved`.

As implemented at S4:

1. The forward sensitivity equations `dS/dt = (∂f/∂x)·S + ∂f/∂θ`, `S(0) = ∂x(0)/∂θ` are generated symbolically and solved with rxode2 (`rtol = 1e-10`, `atol = 1e-12`). Output sensitivities `∂h/∂x·S + ∂h/∂θ` are evaluated on the grid.
2. Known input signals are replaced by a smooth random positive signal (a constant plus three sinusoids with random amplitudes, frequencies and phases), drawn after the Lie points so that the Lie results do not change.
3. Default grid: time 0 and 150 log-spaced times from 0.001 to 200, which covers the rate constants implied by the random parameter values (0.2 to 5).
4. Columns are scaled by the parameter values (log scale); each time point is scaled by its own magnitude (the larger of the output and its sensitivities), and rows that are not finite or below 1e8 × `atol` are dropped (revised at S7 for growing outputs).
5. Rank threshold: relative singular value 1e-7; values between 1e-10 and 1e-7 are ambiguous and give `unresolved`.
6. Reconciliation (`method = "both"`, the default): at every point the two ranks must be equal and the largest principal-angle sine between the two null spaces must be below 1e-4. A parameter keeps its status only if both methods give the same status; otherwise it is `unresolved`. The comparison is returned as the `agreement` table.

**Numerical rules**

- Random points are fixed by a seed; five points by default. A conclusion is drawn only if all points give the same rank.
- Rank is decided from a gap in the singular values; the threshold and the smallest singular values are recorded with the result.
- Parameters are evaluated on a normalized (log) scale.
- If the order of Lie derivatives or the expression size exceeds a limit, the computation stops with `unresolved`. No conclusion is fabricated.

### 5.2 Practical identifiability for a given design

- Expected Fisher information matrix for a typical individual, from sampling times, doses and the residual error model (additive, proportional or combined).
- Eigen-decomposition to report weakly determined directions, the condition number and expected relative standard errors (%).
- Population (between-subject variability) information is considered for v0.2, avoiding overlap with PopED.

As implemented at S5:

1. Output sensitivities come from the same sensitivity equations as S4, solved with rxode2 at the given parameter values; known input signals are given as expressions in `t`.
2. Residual error `sd = sqrt(add^2 + (prop·f)^2)` (`residual_error()`, the nlmixr2 `combined2` form), one per output or shared.
3. Information per observation: `S'S/σ² + 2·prop⁴·f²·S'S/σ⁴` — the second term comes from the dependence of the variance on the parameters. Multiplied by `n` individuals with the same design.
4. Eigen-decomposition on the log-parameter scale; directions with relative eigenvalue below 1e-10 are null. Parameters in a null direction are `non_identifiable` (RSE infinite); otherwise the RSE comes from the (pseudo-)inverse.
5. Status: `identifiable` if RSE ≤ `rse_limit` (default 30%), `unresolved` above it, with a printed warning that the RSEs of all parameters are then approximate and likelihood profiles should be examined.
6. Validation (`validation/fim-simulation.R`): maximum-likelihood estimates from 1,000 simulated data sets per scenario; the fitting routine is cross-checked against derivative-free Nelder–Mead.

### 5.3 Classification of likelihood profiles from actual fits

- For nlmixr2 FOCEi fits, obtain objective function values along a grid of fixed values with `nlmixr2extra::profileFixed()` and classify following Raue et al. (2009):
  - the change in objective function exceeds the threshold (default 3.84) on both sides → `identifiable`;
  - only on one side → `non_identifiable`, with the open direction reported;
  - flat over the range → `non_identifiable`.
- Reconcile with the structural result to distinguish "structurally identifiable but not determined by these data".
- nlmixr2extra is a suggested dependency; without it, only this function is unavailable and says so.

As implemented at S6 (`classify_profiles(fit, which, grid, threshold, steps, flat_tol, structural, map)`):

1. Each side is explored by stepping away from the estimate (default 0.1, 0.25, 0.5, 1, 2 on the estimation scale), stopping once dOFV exceeds the threshold; or along a user grid.
2. Points are evaluated with `nlmixr2extra::profileFixed()`. Because every profile fit is an upper bound of the profile, points more than `flat_tol` (default 1) above the original fit or above the lowest OFV found are re-estimated from a neighbouring profile point (warm start) and the lower value is kept.
3. dOFV is measured from the lowest OFV found; a drop below the original fit is reported.
4. Shapes: bounded, one-sided (open direction), flat (within `flat_tol`), shallow; failed fits give `unresolved`. Approximate interval ends are interpolated on the signed-root scale.
5. With a structural result and a parameter map, each profile is interpreted as "determined by these data", "structurally identifiable but not determined by these data", "structurally non-identifiable", or a conflict.

### 5.4 Which measurement to add

For each candidate supplied by the user (an added output, added time points, an added dose level), report how the structural rank and the expected relative standard errors change. No optimization.

As implemented at S7 (`compare_designs(m, candidates, values, error)`, `experiment(times, events, values, n)`):

1. A candidate is one or more experiments (groups of individuals with the same dosing and sampling). Doses come from the model definition or from an event table (`events`, for example repeated infusions); event doses must be known amounts or rates.
2. The expected information of a candidate is the sum over its experiments; RSEs and statuses as in 5.2.
3. The structural decision is computed once for each distinct set of observed outputs (it does not depend on sampling times or dose levels).
4. Residual error may also be exponential (`residual_error(exp = )`, additive on the log scale, no variance term), matching data fitted as log concentrations.
5. Case studies 2 and 3 are in `case-studies/` (results in `case-studies/README.md`).

## 6. Public functions (v0.1)

| Function | Role |
| --- | --- |
| `pkpd_model(model, outputs, doses, parameters, known)` | Define the model, outputs, doses and known quantities |
| `bolus_dose()`, `infusion_dose()` | Dose specifications (renamed at S7 from `bolus()` and `infusion()`, which clash with rxode2) |
| `structural_identifiability(m, method = c("both", "lie", "sensitivity"), points = 5, seed)` | Structural identifiability |
| `practical_identifiability(m, design, values, error, inputs, n)`, `residual_error(add, prop)` | Practical identifiability from the Fisher information |
| `classify_profiles(fit, which, grid, threshold, structural, map)` | Classification of likelihood profiles (uses nlmixr2extra) |
| `compare_designs(m, candidates, values, error)`, `experiment()` | Comparison of candidate designs (added outputs, sampling times, dose levels) |
| `as_table(x, what)` | Results as tibbles |

Results are S3 objects: `print()` gives a summary and `as_table()` returns machine-readable tables. Each result records the method, random points, thresholds, package version and the statement that the decision is local.

## 7. Validation (core acceptance criteria)

Decisions are checked against **models whose identifiability is known**. The expected results were fixed at milestone S1 from published tables or analytical derivations, before any identifiability code was written. The full record, with sources, equations and derivations, is [validation/reference-models.md](validation/reference-models.md).

**Acceptance set (pass/fail)**

| ID | Model | Expected result | Source |
| --- | --- | --- | --- |
| A1 | One-compartment IV bolus | CL and V identifiable | Derivation |
| A2 | One-compartment oral, F unknown | F, V, CL not identifiable alone; `V/F`, `CL/F` identifiable; ka locally identifiable | Janzén et al. 2016 (citing Cheung et al. 2013); derivation |
| A3 | One-compartment oral, F = 1 | Locally identifiable; two global solutions (flip-flop) | Derivation |
| A4 | Two-compartment IV, central observed | CL, V1, Q, V2 identifiable | Derivation |
| A5 | Michaelis–Menten elimination, IV | Vmax, Km, V identifiable | Derivation |
| A6 | Parallel linear and Michaelis–Menten elimination, IV | CL, Vmax, Km, V identifiable | Derivation |
| A7 | Direct effect, steady-state binding, linear (known `Cp`) | `Rtot·ke`, Kd identifiable; Rtot, ke not alone | Janzén et al. 2016, Model 1 |
| A8 | Indirect stimulation, steady-state binding (known `Cp`) | `Rtot·ke`, kin, kout, Kd identifiable | Janzén et al. 2016, Model 3 |
| A9 | Effect compartment, steady-state binding, linear (known `Cp`) | `Rtot·ke`, ke0, Kd identifiable | Janzén et al. 2016, Model 9 |
| A10 | Effect compartment, dynamic binding, linear (known `Cp`) | `Rtot·ke`, ke0, kon, koff identifiable; all identifiable with Rtot fixed | Janzén et al. 2016, Model 13 |

**Exploratory set (reported, not pass/fail):** full target-mediated drug disposition and its approximations (Eudy et al. 2015; conditions not verifiable from the abstract), and Janzén et al. 2016 Models 2, 6, 10 and 14 (inconsistent table entries). The draft expectation that full TMDD parameters require a target measurement was withdrawn at S1.

In addition:

- Agreement between the primary (Lie) and secondary (sensitivity) methods is reported for all models.
- Practical identifiability: on synthetic data with known values, expected relative standard errors from the Fisher information agree with the spread of repeated estimates.
- Profile classification: on nlmixr2 fits to synthetic data constructed to be identifiable or not (for example the oral model with F unknown), the classification is correct.

## 8. Case studies (vignettes)

Using only public data and connecting to published analyses:

1. Decisions for the ten acceptance models above.
2. Nimotuzumab target-mediated disposition: which parameters are determined by the observations in the public `nimoData` dataset (nlmixr2data), and which measurement would determine the rest.
3. PSA decline and regrowth: how the identifiability of the regrowth rate depends on the length of observation.

## 9. Quality and operation

- R package with testthat tests and GitHub Actions R CMD check on Windows and Linux.
- Runs fully offline with no external communication. Random numbers are seeded and results are reproducible.
- License GPL-3 (rxode2 and symengine are GPL).
- All documentation (README, vignettes, help pages, NEWS) is in English.
- Commits use the noreply identity. Before any publication, the private-term check is run with zero hits. Only public or synthetic data are used.
- Errors and warnings carry codes (for example `PKI001 MODEL_PARSE_FAILED`, `PKI010 SENSITIVITY_SOLVE_FAILED`).
- The repository is made public after milestone S3.

## 10. Milestones

| Stage | Content | Done when |
| --- | --- | --- |
| S1 | Fix expected results for the reference models from primary literature | **Done 2026-10-05**: `validation/reference-models.md` (10 acceptance models, exploratory set) |
| S2 | `pkpd_model()`: extract equations, outputs and inputs (including known input signals) from rxode2 models | **Done 2026-10-05**: all acceptance models (and the Rtot-known variant) load; 9 tests, 51 expectations; R CMD check OK |
| S3 | Structural identifiability (Lie method) | **Done 2026-10-05**: all 10 acceptance models (and the Rtot-known variant) agree with expectations at 4 seeds; 17 tests, 196 expectations; R CMD check OK. Repository public 2026-10-05 |
| S4 | Sensitivity check and reconciliation of the two methods | **Done 2026-10-05**: the two methods agree at every point for all acceptance models (and the Rtot-known variant) at 4 seeds, and for exploratory TMDD; the sensitivity method alone reproduces the expected results; disagreement gives `unresolved` (tested); 19 tests, 254 expectations; R CMD check 0 errors, 0 warnings, 0 notes |
| S5 | Practical identifiability (Fisher information) | **Done 2026-10-05**: information matches the analytic result (relative difference 1e-11); in four acceptance scenarios (A1, A3, A4, A10 with Rtot known) expected RSEs agree with the spread of 1,000 repeated estimates within 0.96–1.08; limitations recorded (RSE 30–40%: skewed estimates, spread underestimated) and the default RSE limit set to 30%; 27 tests, 280 expectations; R CMD check 0 errors, 0 warnings, 0 notes |
| S6 | Profile classification (nlmixr2extra) | **Done 2026-10-05**: on three synthetic FOCEi fits (identifiable; F unknown; Emax with concentrations far below EC50) all ten parameters are classified as expected, including shape and open direction, and profile intervals cover the true values; joined with structural results; 35 tests, 309 expectations; R CMD check 0 errors, 0 warnings, 0 notes |
| S7 | Comparison of added measurements | **Done 2026-10-05**: case study 2 (nimoData, full TMDD: cl, v, koff determined; kon, kdeg, r0, kint flat; total target would determine 6 of 7) and case study 3 (PSA: regrowth rate determined only after about 12 months) produce results; event-table dosing and exponential error checked against analytic information; sensitivity check revised for growing outputs (S4 rerun, unchanged); 42 tests; R CMD check 0 errors, 0 warnings, 0 notes |
| S8 | Vignettes, CI, README, version 0.1.0 | **Done 2026-10-05**: three vignettes (reference models, nimotuzumab TMDD, PSA regrowth; heavy fits read from `inst/extdata`); GitHub Actions R CMD check on Linux and Windows; README and NEWS; version 0.1.0; R CMD check with vignettes 0 errors, 0 warnings, 0 notes; both case-study scripts rerun and every result table identical byte for byte |

## 11. Main risks

| Risk | Mitigation |
| --- | --- |
| Expression swell in symbolic computation (for example TMDD) | Limits on order and expression size, `unresolved` beyond them; numerical sensitivity method alongside |
| Wrong rank decisions at random points | Several points, normalization, recorded singular values; `unresolved` when the methods disagree |
| Confusing local and global identifiability | "Local" stated in every result; flip-flop included among reference models |
| Changes in rxode2 internals | Minimal use of internal functions, versions recorded, detected by tests |
| Scope creep | v0.1 fixed to section 2; variability, population information and global decisions deferred to v0.2 or later |

## 12. Decisions taken

- Name: `pkident` (PK identifiability; not used on CRAN, Bioconductor, GitHub or PyPI as of 2026-10-05). Earlier working names: pmxident, pkpdIdentifiability.
- Public repository after S3.
- All documentation in English.
- No comparison with Julia tools in validation.
