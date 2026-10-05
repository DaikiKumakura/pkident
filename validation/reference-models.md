# Reference models and expected results (milestone S1)

Fixed on 2026-10-05, before any identifiability code was written. Each expected result comes from a published table or an analytical derivation given here. Results that could not be verified are kept out of the acceptance set.

Notation: states are amounts unless stated; `D` is a known dose; outputs are observed continuously and without noise (the definition of structural identifiability). "Local" means identifiable in a neighbourhood of the true values; "global" means unique over the parameter space. pkident decides local identifiability, so global results below are compared at the local level (global identifiability implies local identifiability; parameters that enter only through a product are locally unidentifiable as well).

## Sources

| Key | Reference | Access used |
| --- | --- | --- |
| Janzen2016 | Janzén DLI, Bergenholm L, Jirstrand M, Parkinson J, Yates J, Evans ND, Chappell MJ. Parameter identifiability of fundamental pharmacodynamic models. *Front Physiol*. 2016;7:590. doi:10.3389/fphys.2016.00590 (PMCID PMC5136565) | Full text (open access): Tables 1–2, equations 1, 6–8 |
| Cheung2013 | Cheung SYA, Yates JWT, Aarons L. The design and analysis of parallel experiments to produce structurally identifiable models. *J Pharmacokinet Pharmacodyn*. 2013;40:93–100. doi:10.1007/s10928-012-9291-z | Abstract; cited by Janzen2016 for the oral model |
| Eudy2015 | Eudy RJ, Riggs MM, Gastonguay MR. A priori identifiability of target-mediated drug disposition models and approximations. *AAPS J*. 2015;17:1280–1284. doi:10.1208/s12248-015-9795-8 | Abstract only (full text not open access) |
| Godfrey1994 | Godfrey KR, Chapman MJ, Vajda S. Identifiability and indistinguishability of nonlinear pharmacokinetic models. *J Pharmacokinet Biopharm*. 1994;22:229–251. doi:10.1007/BF02353330 | Abstract only; context for Michaelis–Menten models |

## Acceptance set

### A1. One-compartment IV bolus

`dA/dt = −(CL/V)·A`, `A(0) = D`, `y = A/V`.

Expected: **CL and V identifiable** (global).
Derivation: `y(t) = (D/V)·exp(−(CL/V)t)`. `y(0) = D/V` gives V; the decay rate gives CL/V and hence CL.

### A2. One-compartment oral, bioavailability unknown

`dAd/dt = −ka·Ad`, `Ad(0) = F·D`; `dAc/dt = ka·Ad − (CL/V)·Ac`, `Ac(0) = 0`; `y = Ac/V`. Parameters: ka, CL, V, F.

Expected: **F, V and CL not identifiable individually; `V/F` and `CL/F` identifiable combinations; ka locally identifiable** (globally two solutions; see A3).
Source: Janzen2016 (equation 1 and the statement that only the fraction F/V can be identified, citing Cheung2013). Derivation: `y(t) = (F·D/V)·ka/(ka − k)·(exp(−kt) − exp(−ka·t))` with `k = CL/V`; F and V enter only through `F/V`, and CL only through `k = (CL/F)/(V/F)`.

### A3. One-compartment oral, F = 1 known

As A2 with F fixed to 1. Parameters: ka, CL, V.

Expected: **ka, CL and V locally identifiable; not globally identifiable (flip-flop: two solutions).**
Derivation: with `k = CL/V`, `y(t) = (D/V)·ka/(ka − k)·(exp(−kt) − exp(−ka·t))`. The set `(ka', k', V') = (k, ka, V·k/ka)` gives the same `y(t)` for all t, and no other set does (the two exponents are determined as an unordered pair, and the amplitude then fixes V). The solutions are isolated, so the parameters are locally but not globally identifiable. The test checks the local decision and that every result states that it is local.

### A4. Two-compartment IV bolus, central compartment observed

`dA1/dt = −(CL/V1 + Q/V1)·A1 + (Q/V2)·A2`, `A1(0) = D`; `dA2/dt = (Q/V1)·A1 − (Q/V2)·A2`, `A2(0) = 0`; `y = A1/V1`. Parameters: CL, V1, Q, V2.

Expected: **CL, V1, Q and V2 identifiable** (global).
Derivation: `y(t) = P·exp(−αt) + R·exp(−βt)` determines `P, R, α, β`. Then `V1 = D/(P + R)`; `k21 = (Pβ + Rα)/(P + R)`; `k10 = αβ/k21`; `k12 = α + β − k10 − k21`; hence `CL = k10·V1`, `Q = k12·V1`, `V2 = Q/k21`. The map is one-to-one.

### A5. One-compartment Michaelis–Menten elimination, IV bolus

`dA/dt = −Vmax·(A/V)/(Km + A/V)`, `A(0) = D`, `y = A/V`. Parameters: Vmax, Km, V.

Expected: **Vmax, Km and V identifiable.**
Derivation: `y(0) = D/V` gives V. Then `dy/dt = −(Vmax/V)·y/(Km + y)` is observed as a function of y over a continuous range of concentrations, which determines `Vmax/V` and `Km`, hence Vmax. (Godfrey1994 analyses two-compartment Michaelis–Menten variants; used as context only.)

### A6. One-compartment parallel linear and Michaelis–Menten elimination, IV bolus

`dA/dt = −(CL/V)·A − Vmax·(A/V)/(Km + A/V)`, `A(0) = D`, `y = A/V`. Parameters: CL, Vmax, Km, V.

Expected: **CL, Vmax, Km and V identifiable.**
Derivation: V from `y(0)`. `dy/dt = −k·y − a·y/(Km + y)` with `k = CL/V`, `a = Vmax/V`. On any interval of concentrations, the functions `y` and `y/(Km + y)` are linearly independent for `Km > 0`, so `k`, `a` and `Km` are determined. Practical identifiability requires concentrations spanning the nonlinear range; that is not part of this structural test.

### A7–A10. Pharmacodynamic models with a known plasma concentration input (Janzen2016)

Plasma concentration `Cp(t)` is a **known input signal**. These four models were chosen because their entries in Janzen2016 Table 2 are internally consistent. Results are for the fixed-effects versions with Rtot (total receptors) estimated.

| ID | Janzen2016 model | Equations | Expected (Table 2) |
| --- | --- | --- | --- |
| A7 | Model 1: direct, steady-state binding, linear transduction | `E = ke·Rtot·Cp/(Kd + Cp)` | **Rtot and ke not identifiable individually; `Rtot·ke` and Kd identifiable** |
| A8 | Model 3: direct, steady-state binding, indirect (stimulation of kin) | `dE/dt = kin·(1 + ke·Rtot·Cp/(Kd + Cp)) − kout·E`, `E(0) = kin/kout` | **Rtot, ke not identifiable individually; `Rtot·ke`, kin, kout, Kd identifiable** |
| A9 | Model 9: effect compartment, steady-state binding, linear | `dCe/dt = ke0·(Cp − Ce)`, `Ce(0) = 0`; `E = ke·Rtot·Ce/(Kd + Ce)` | **Rtot, ke not identifiable individually; `Rtot·ke`, ke0, Kd identifiable** |
| A10 | Model 13: effect compartment, dynamic binding, linear | `dCe/dt = ke0·(Cp − Ce)`; `dRC/dt = kon·(Rtot − RC)·Ce − koff·RC`, `RC(0) = 0`; `E = ke·RC` | **Rtot, ke not identifiable individually; `Rtot·ke`, ke0, kon, koff identifiable.** With Rtot fixed (the paper's worked example, equations 6–8): ke0, ke, kon, koff globally identifiable |

Notes on A8: Janzen2016 Table 1 prints the initial condition as `E(0) = kout/kin`. The steady-state baseline of this model is `kin/kout`, so the printed form is treated as a typographical inversion and `E(0) = kin/kout` is used. The expected result does not depend on this choice of form as long as the baseline is at steady state.

Note on A10: rescaling `r = RC/Rtot` gives `dr/dt = kon·(1 − r)·Ce − koff·r` and `E = (ke·Rtot)·r`, which confirms that Rtot and ke enter only through their product while kon and koff remain identifiable.

## Exploratory set (reported, not pass/fail)

| ID | Model | Status of the reference | Use |
| --- | --- | --- | --- |
| E1 | Full target-mediated drug disposition (and QE/RB, QSS, MM approximations) | Eudy2015 abstract states the full model and approximations are a priori identifiable "regardless of whether observations were taken from a single or multiple compartments". The observed quantity and the assumptions about known initial conditions could not be verified without the full text | Case study; results reported with their own assumptions, not used as acceptance |
| E2 | Janzen2016 Models 2, 6, 10, 14 (sigmoid transduction) | Table 2 entries list parameters not present in the corresponding equations (for example ke for Model 2, kon/koff for Model 10, kin for Model 14) | Excluded from acceptance until the inconsistency is resolved (for example from the supplementary material) |

## Correction to the specification draft

The draft expected that full TMDD parameters are not determined without a target measurement. Eudy2015 reports a priori identifiability from a single compartment. TMDD therefore moves to the exploratory set, and the draft expectation is withdrawn.

## Implication for the specification

The PD reference models need a **known input signal** (`Cp(t)`), not only bolus, infusion and depot doses. pkident v0.1 therefore supports known inputs whose time derivatives are treated as known quantities in the Lie-derivative method and as a known forcing function in the sensitivity method.

## Results at milestone S3 (2026-10-05)

`structural_identifiability()` with default settings (5 random points, seed 1). Every acceptance model agrees with its expected result; the decisions are the same with seeds 2, 3 and 4.

| ID | Decision | Derivative order | Rank | Identifiable | Not identifiable alone | Identifiable combinations reported |
| --- | --- | ---: | ---: | --- | --- | --- |
| A1 | agrees | 1 | 2/2 | CL, V | — | — |
| A2 | agrees | 8 | 3/4 | ka | CL, F, V | CL/F, CL/V (spans V/F) |
| A3 | agrees | 3 | 3/3 | CL, ka, V | — | — |
| A4 | agrees | 3 | 4/4 | CL, Q, V1, V2 | — | — |
| A5 | agrees | 2 | 3/3 | KM, V, VMAX | — | — |
| A6 | agrees | 3 | 4/4 | CL, KM, V, VMAX | — | — |
| A7 | agrees | 5 | 2/3 | Kd | Rtot, ke | Rtot*ke |
| A8 | agrees | 8 | 4/5 | Kd, kin, kout | Rtot, ke | Rtot*ke |
| A9 | agrees | 7 | 3/4 | Kd, ke0 | Rtot, ke | Rtot*ke |
| A10 | agrees | 9 | 4/5 | ke0, koff, kon | Rtot, ke | Rtot*ke |
| A10, Rtot known | agrees | 5 | 4/4 | ke, ke0, koff, kon | — | — |

For A2 the reported basis (CL/F, CL/V) differs from the form in the expected result (V/F, CL/F) but spans the same set of combinations; the test checks that V/F is orthogonal to the null space.

Exploratory E1 (full TMDD; free drug observed; bolus `L(0) = DOSE/V`; target at its steady-state baseline `R(0) = ksyn/kdeg`; complex initially absent): all seven parameters (V, kel, kon, koff, ksyn, kdeg, kint) locally identifiable at derivative order 9. This agrees with the abstract of Eudy2015; the conditions of that paper could not be verified, so the case remains exploratory.

## Results at milestone S4 (2026-10-05)

`structural_identifiability(method = "both")`: the Lie method is checked at each of the 5 random points by the rank and null space of the output sensitivity matrix (rxode2, 151 time points from 0 to 200; known input signals replaced by a random smooth signal).

| ID | Agreement (points, seeds 1–4) | Largest null-space distance | Smallest retained sensitivity ratio | Largest discarded ratio |
| --- | --- | ---: | ---: | ---: |
| A1 | 20/20 | 0 | 1.8e-1 | — |
| A2 | 20/20 | 0 | 1.2e-2 | 6.7e-15 |
| A3 | 20/20 | 0 | 1.2e-2 | — |
| A4 | 20/20 | 0 | 3.0e-4 | — |
| A5 | 20/20 | 0 | 2.9e-3 | — |
| A6 | 20/20 | 0 | 1.6e-6 | — |
| A7 | 20/20 | 0 | 6.8e-3 | 1.1e-16 |
| A8 | 20/20 | 0 | 4.6e-4 | 1.3e-16 |
| A9 | 20/20 | 0 | 7.4e-3 | 1.9e-16 |
| A10 | 20/20 | 1.5e-8 | 2.0e-3 | 3.7e-16 |
| A10, Rtot known | 20/20 | 0 | 3.7e-3 | — |

Every decision is the same as at S3, and the sensitivity method alone also reproduces every expected result. Ratios are singular values relative to the largest, on the log-parameter scale; the threshold is 1e-7 with an ambiguous band down to 1e-10. The structural null directions lie at about 1e-16, at least nine orders of magnitude below the smallest retained value.

The smallest retained values for A6 (1.6e-6, parallel linear and Michaelis–Menten elimination) and for exploratory E1 (1.6e-7, full TMDD) are close to the threshold. These directions are weakly determined at some random points (practical, not structural, weakness), and they are the reason for the ambiguous band: a value falling inside it gives `unresolved` rather than a decision.

Exploratory E1 (full TMDD, conditions as at S3): both methods give rank 7 of 7 at all 5 points, so all seven parameters are locally identifiable.

A test with only two time points (A4, `times = c(0, 0.001)`) checks the reconciliation rule: the sensitivity rank is too low, the methods disagree, and the result is `unresolved` with no parameter reported as non-identifiable.
