# pkident

Structural and practical identifiability for pharmacokinetic and pharmacodynamic ODE models.

**Status: in development (not released).** Milestones S1–S3 are complete: models are read from rxode2, and local structural identifiability is decided for all reference models.

```r
library(pkident)

m <- pkpd_model(
  "d/dt(depot)   = -ka*depot
   d/dt(central) = ka*depot - CL/V*central
   cp            = central/V",
  outputs = "cp",
  doses   = bolus("depot", "F*DOSE"),
  known   = "DOSE"
)

structural_identifiability(m)
#> <pkident structural identifiability>
#> Decision is LOCAL (Lie derivatives at the known initial state; 5 random points, derivative order up to 8).
#>   CL           combination_only
#>   F            combination_only
#>   V            combination_only
#>   ka           identifiable
#> Identifiable combinations:
#>   CL/F
#>   CL/V
```

Decisions are local: a parameter reported as identifiable may still have a finite number of alternative values elsewhere in parameter space (for example the flip-flop of absorption and elimination).

See `SPEC.md` for the specification and `validation/reference-models.md` for the reference models, their expected results and the results obtained.

## License

GPL-3
