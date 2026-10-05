# Practical identifiability: validation (milestone S5)

Run on 2026-10-05 with `Rscript validation/fim-simulation.R` (pkident 0.0.0.9000, rxode2 5.1.8, symengine 0.2.13). Raw results: `fim-simulation-results.csv`.

## Checks

1. **Analytic check.** For the one-compartment IV bolus model, the Fisher information from `practical_identifiability()` equals the closed-form information (including the term from a proportional error variance) with relative difference about 1e-11, for combined, proportional-only and additive-only error (`tests/testthat/test-practical.R`).
2. **Repeated estimation.** For each scenario, 1,000 data sets are simulated at the true values with the residual error model and fitted by maximum likelihood (residual error known; Fisher scoring with Levenberg–Marquardt damping on log parameters). The expected RSE is compared with the standard deviation of the log estimates (empirical RSE) and with 100 × the median absolute deviation (robust RSE, insensitive to skewed tails).
3. **Fitting routine.** Five A3 data sets refitted with derivative-free Nelder–Mead give the same estimates (largest relative difference 1.1e-7).

Monte Carlo uncertainty of an empirical RSE from 1,000 data sets is about 2%.

## Acceptance scenarios (all expected RSEs below 15%)

| Scenario | Design and error | Parameter | Expected RSE % | Empirical RSE % | Ratio |
| --- | --- | --- | ---: | ---: | ---: |
| A1 one-compartment IV bolus | 8 samples 0.25–24 h; add 0.05, prop 10% | CL | 4.69 | 4.86 | 1.04 |
| | | V | 4.75 | 4.55 | 0.96 |
| A3 one-compartment oral, F known | same | CL | 4.97 | 5.17 | 1.04 |
| | | V | 8.04 | 8.34 | 1.04 |
| | | ka | 13.44 | 14.27 | 1.06 |
| A4 two-compartment IV bolus | 10 samples 0.1–48 h; add 0.02, prop 10% | CL | 4.82 | 4.88 | 1.01 |
| | | Q | 11.54 | 11.72 | 1.02 |
| | | V1 | 6.14 | 6.21 | 1.01 |
| | | V2 | 11.61 | 12.05 | 1.04 |
| A10 effect compartment, dynamic binding, Rtot known; `Cp = 10·exp(−0.2t)` | hourly 0–24 h; add 0.02, prop 1% | ke | 1.17 | 1.26 | 1.08 |
| | | ke0 | 6.50 | 6.56 | 1.01 |
| | | koff | 8.30 | 8.44 | 1.02 |
| | | kon | 6.50 | 6.50 | 1.00 |

All 4,000 fits converged. Every ratio is between 0.96 and 1.08: the expected RSEs agree with the spread of repeated estimates.

## Limitations found (reported, not pass/fail)

| Scenario | Parameter | Status | Expected RSE % | Empirical RSE % | Robust RSE % |
| --- | --- | --- | ---: | ---: | ---: |
| A10 as above with add 0.1, prop 5% (952 of 1,000 fits converged) | ke | identifiable | 5.83 | 18.8 | 7.55 |
| | ke0 | unresolved | 32.5 | 44.5 | 32.0 |
| | koff | unresolved | 41.4 | 53.6 | 38.6 |
| | kon | unresolved | 32.4 | 42.3 | 32.1 |
| A5 Michaelis–Menten, concentrations mostly below KM (797 of 1,000 converged) | KM | unresolved | 150 | 142 | 121 |
| | V | identifiable | 5.44 | 4.75 | 4.54 |
| | VMAX | unresolved | 97.6 | 90.8 | 69.7 |

With RSEs of 30–40%, the centre of the distribution of estimates still matches the expected RSE (robust RSE), but the estimates have skewed tails and the standard deviation is 30% larger. The tails also widen the spread of `ke`, whose own expected RSE is small (5.8% expected, 18.8% empirical). The expected information is a local, linear approximation, and this is where it stops being reliable.

Consequences in the implementation:

- The default `rse_limit` is 30%: parameters above it are `unresolved`, not `identifiable`.
- When any parameter is `unresolved`, the printed result states that the RSEs of all parameters may be underestimated and that likelihood profiles should be examined (milestone S6).
