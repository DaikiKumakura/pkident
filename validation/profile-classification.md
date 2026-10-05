# Profile classification: validation (milestone S6)

Run on 2026-10-05 with `Rscript validation/profile-classification.R` (nlmixr2est 7.1.0, nlmixr2extra 5.2.1, rxode2 5.1.8; FOCEi). Raw results: `profile-classification-results.csv` (classification) and `profile-classification-profiles.csv` (every profile point).

## Synthetic cases

Data are simulated so that the answer is known by construction (`tests/testthat/helper-profile-cases.R`; 24 subjects each).

| Case | Data | Fitted model | Expected |
| --- | --- | --- | --- |
| oral | One-compartment oral, dose 100, 9 samples 0.25–24 h, proportional error 10%, between-subject variability on ka and CL | ka, CL, V (F = 1) | All identifiable |
| oral_f | Same data | ka, CL, V and F estimated | ka identifiable; CL, V, F flat (only CL/F and V/F are determined) |
| emax_low | Direct Emax effect with known concentrations 0.02–0.5, far below EC50 = 20; additive error 0.3 | E0, Emax, EC50 | E0 identifiable; Emax and EC50 bounded below, open upwards |

## Results

| Case | Parameter | Expected | Result | Shape | Approx. interval (log scale) | Structural decision | Interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| oral | tka | identifiable | identifiable | bounded | −0.103 to 0.156 (truth 0) | identifiable | determined by these data |
| | tcl | identifiable | identifiable | bounded | 0.619 to 0.797 (truth 0.693) | identifiable | determined by these data |
| | tv | identifiable | identifiable | bounded | 2.947 to 3.004 (truth 2.996) | identifiable | determined by these data |
| oral_f | tka | identifiable | identifiable | bounded | −0.121 to 0.175 | identifiable | determined by these data |
| | tcl | non-identifiable, flat | non_identifiable | flat | — | combination_only | structurally non-identifiable |
| | tv | non-identifiable, flat | non_identifiable | flat | — | combination_only | structurally non-identifiable |
| | tf | non-identifiable, flat | non_identifiable | flat | — | combination_only | structurally non-identifiable |
| emax_low | te0 | identifiable | identifiable | bounded | 2.254 to 2.416 | identifiable | determined by these data |
| | temax | open upwards | non_identifiable | one-sided, open upper | lower end 0.683 | identifiable | structurally identifiable but not determined by these data |
| | tec50 | open upwards | non_identifiable | one-sided, open upper | lower end −0.692 | identifiable | structurally identifiable but not determined by these data |

All ten classifications agree with the expected results, and the three profile intervals of the identifiable case cover the true values. The structural decisions come from `structural_identifiability()` on A3, A2 and `effect = E0 + Emax*Cp/(EC50 + Cp)` with `Cp` a known input; joining them separates structural non-identifiability (oral_f) from data that do not determine structurally identifiable parameters (emax_low).

Explored ranges: ±2 on the log scale where the threshold was not reached (a factor of 7.4 in each direction). Statements about open directions apply to this range only.

## Problems found during validation, and the changes they led to

1. **Spurious high values along flat directions.** With every profile fit started from the original estimates (`nlmixr2extra::profileFixed()`), two of thirty points along the flat directions of oral_f gave dOFV ≈ 155: far from the estimate, the optimiser did not reach the ridge. Every profile fit is an upper bound of the profile, so failures can only overstate dOFV. Points more than `flat_tol` above the original fit or above the lowest OFV found are now re-estimated from the estimates of a neighbouring profile point (warm start), keeping the lower value. Nine points were improved in this way in the final run.
2. **Original fit above the minimum.** Along the flat directions of oral_f, profiling found OFVs 0.92–0.95 below the original fit. dOFV is now measured from the lowest value found and the drop is reported (`ofv_drop`); previously such profiles were marked `unresolved`.
3. **Interval ends.** Linear interpolation of dOFV between the estimate and a point with dOFV 47 gave an interval for tv (2.97 to 2.98) that missed the true value. Interpolation is now on the signed-root scale, which is exact for a quadratic profile; the interval (2.947 to 3.004) covers the true value. The intervals remain approximate; `profile(fit, method = "llp")` gives precise ones.
4. **Case design.** The first Emax design (concentrations up to 4, EC50 = 20) bounded Emax on both sides within the explored range; the data were more informative than intended. Concentrations were lowered to at most 0.5.

## Limits

- The flat/shallow distinction depends on `flat_tol` (default 1). Repeated FOCEi estimation along the flat ridge varied by up to 0.95 here, close to the tolerance; the status (`non_identifiable`) does not depend on it.
- Run time: about 5 minutes for four parameters with flat profiles (10 profile fits each, plus warm starts).
