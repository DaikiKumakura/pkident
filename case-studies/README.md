# Case studies (milestone S7)

Run on 2026-10-05 from the package root (pkident 0.0.0.9000, rxode2 5.1.8, nlmixr2est 7.1.0, nlmixr2extra 5.2.1). Result tables are in `inst/extdata/` (installed with the package and used by the vignettes). The vignettes present the same analyses.

## Case study 2: nimotuzumab target-mediated disposition (`nimotuzumab-tmdd.R`)

**Question.** Which parameters of a full target-mediated drug disposition (TMDD) model are determined by the public `nimoData` observations (nlmixr2data; source: Rodriguez-Vera et al. 2015), and which added measurement would determine the rest?

**Data and model.** 12 patients, 50–400 mg, ten repeated infusions each, 321 free-drug observations (DV is the log concentration; its unit is not documented, so the volume is in dose units per concentration unit). Full TMDD model (drug amount `A`, free target `R`, complex `P`; target at its baseline `r0` before the first dose; parameters cl, v, kon, koff, kdeg, r0, kint), fitted with nlmixr2 FOCEi, between-subject variability on cl and v, additive error on the log scale (SD 0.668).

**1. Structural.** With free drug observed alone, all seven parameters are locally identifiable (both methods agree at all points). Structure is not the obstacle.

**2. What these data determine (likelihood profiles, `classify_profiles()`).**

| Parameter | Estimate (log) | Profile | Interpretation |
| --- | ---: | --- | --- |
| cl | −5.44 | bounded, −6.91 to −4.88 | determined by these data |
| v | 0.30 | bounded, −0.04 to 0.66 | determined by these data |
| koff | −8.36 | bounded, −9.08 to −7.26 | determined by these data |
| kon | −6.07 | flat over ±2 | structurally identifiable but not determined by these data |
| kdeg | −1.29 | flat over ±2 | structurally identifiable but not determined by these data |
| r0 | 1.01 | flat over −2 to +2 | structurally identifiable but not determined by these data |
| kint | −12.25 | flat over ±2 | structurally identifiable but not determined by these data |

The four flat parameters stay within dOFV 0.9 of the minimum over a factor of 7.4 in each direction. nlmixr2's covariance step reports standard errors of 0.17–0.48 on the log scale for these parameters; the profiles show that these standard errors do not describe the uncertainty.

**3. Which measurement would determine the rest (`compare_designs()`).** Expected RSE (%) at the fitted values (fixed effects of a typical patient; the actual infusions and sampling times of each of the 12 patients). Target measurements are assumed to have an error of 0.2 on the log scale.

| Parameter | Current design | + total target | + free target | + complex | + 10 mg group (3 patients) |
| --- | ---: | ---: | ---: | ---: | ---: |
| cl | 35 | 7.9 | 9.5 | 7.6 | 29 |
| v | 9.0 | 5.4 | 5.9 | 5.4 | 7.7 |
| kon | 2,838 | 9.4 | 13 | 79 | 2,574 |
| koff | 242 | 23 | 89 | 22 | 173 |
| kdeg | 2,849 | 7.6 | 16 | 78 | 2,584 |
| r0 | 2,830 | 5.7 | 1.9 | 77 | 2,568 |
| kint | 44,217 | 1,055 | 16,912 | 1,056 | 31,796 |
| RSE ≤ 30% | 1/7 | 6/7 | 5/7 | 3/7 | 2/7 |

**Answer.** These data determine clearance, volume and (loosely) koff. Adding a lower dose level does not help. Measuring total target at the same times would determine kon, kdeg and r0 as well (6 of 7 parameters with RSE ≤ 30%); free target determines five. No candidate determines kint: at the estimate (4.8e-6 per hour) internalisation of the complex has almost no effect on any output, so its value is not identifiable in practice whatever is measured; a design aimed at kint would need a different design point or a different model.

**Agreement between the methods.** The profiles and the expected RSEs agree on which parameters are determined (cl, v and koff bounded; kon, kdeg, r0 and kint flat versus RSEs above 2,800%). They differ in degree: the expected RSE of koff (242%) exceeds what its bounded profile suggests, and that of cl (35%) is smaller. The expected information is for fixed effects of a typical patient and ignores between-subject variability, and it is a local approximation; the profiles come from the actual mixed-effects fit.

## Case study 3: PSA decline and regrowth (`psa-regrowth.R`)

**Question.** How does the identifiability of the regrowth rate depend on the length of observation?

**Model and values.** `S' = −d·S`, `R' = g·R`, `PSA = S + R`, the ODE form of the biexponential decline-and-regrowth model used for PSA and tumour size (for example Stein et al. 2008). The values are **illustrative**, not estimates: PSA 20 ng/mL at the start (s0 = 19.5, r0 = 0.5), d = 0.8 per month, g = 0.15 per month; nadir at 5.6 months. Error 15% proportional plus 0.05 ng/mL; one patient sampled monthly.

**Structural.** All four parameters are locally identifiable from PSA alone (both methods agree). This case exposed a weakness in the sensitivity check for growing outputs, which was corrected (see `validation/reference-models.md`).

**Expected RSE (%) by length of observation.**

| Parameter | 3 months | 6 months | 9 months | 12 months | 18 months | 24 months |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| d (decline) | 118 | 20 | 13 | 10 | 8.9 | 8.4 |
| g (regrowth) | 5,298 | 200 | 47 | 21 | 7.9 | 4.4 |
| r0 | 3,178 | 183 | 56 | 30 | 15 | 11 |
| s0 | 78 | 13 | 13 | 13 | 12 | 12 |

**Check by repeated estimation** (1,000 simulated patients each): at 12 months the ratio of empirical to expected RSE is 1.04–1.06 for all parameters; at 24 months 1.01–1.06.

**Answer.** The regrowth rate is not determined until PSA has been followed well past the nadir: with these values its RSE falls below 30% only at about 12 months (twice the time to nadir), while the decline rate is determined by 6 months. An analysis that stops near the nadir can describe the decline but not the regrowth.
