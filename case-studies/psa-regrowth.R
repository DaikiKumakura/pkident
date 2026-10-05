# Case study 3 (milestone S7): PSA decline and regrowth.
#
# How does the identifiability of the regrowth rate depend on the length of
# observation? Model: a treatment-sensitive PSA component that declines and
# a resistant component that grows,
#   S' = -d S, R' = g R, PSA = S + R, S(0) = s0, R(0) = r0,
# the ODE form of the biexponential decline-and-regrowth model used for PSA
# and tumour size (for example Stein et al., Clin Cancer Res 2008).
#
# Parameter values are ILLUSTRATIVE, not estimates from data: PSA 20 ng/mL at
# the start of treatment, 2.5% from the resistant component, sensitive
# half-life about 3.5 weeks (d = 0.8 per month), resistant doubling time about
# 4.6 months (g = 0.15 per month). Residual error 15% proportional plus
# 0.05 ng/mL additive (assumed). One patient sampled monthly from month 0.
#
# Run from the package root: Rscript case-studies/psa-regrowth.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-fim-simulation.R")
out_dir <- "case-studies/results"
dir.create(out_dir, showWarnings = FALSE)

m <- pkpd_model(
  "d/dt(S) = -d*S
   d/dt(R) = g*R
   psa = S + R",
  outputs = "psa", initial = list(S = "s0", R = "r0")
)
values <- c(s0 = 19.5, r0 = 0.5, d = 0.8, g = 0.15)
error <- residual_error(add = 0.05, prop = 0.15)
nadir <- log(values[["d"]] * values[["s0"]] / (values[["g"]] * values[["r0"]])) / (values[["d"]] + values[["g"]])
cat(sprintf("PSA nadir at %.1f months\n", nadir))

print(structural_identifiability(m))

lengths <- c(3, 6, 9, 12, 18, 24)
candidates <- stats::setNames(lapply(lengths, function(L) experiment(seq(0, L, 1))),
                              paste0(lengths, " months"))
cmp <- compare_designs(m, candidates, values = values, error = error, structural = FALSE)
print(cmp)
utils::write.csv(cmp$summary, file.path(out_dir, "psa-designs-summary.csv"), row.names = FALSE)
utils::write.csv(cmp$parameters, file.path(out_dir, "psa-designs-parameters.csv"), row.names = FALSE)

# Check of the expected RSEs by repeated estimation for two lengths
check <- list()
for (L in c(12, 24)) {
  tt <- seq(0, L, 1)
  p <- practical_identifiability(m, tt, values, error)
  sim <- simulate_and_fit(m, tt, values, error, reps = 1000, seed = 2026)
  cmp_l <- compare_rse(p, sim)
  le <- log(sim$estimates[sim$converged, , drop = FALSE])
  cmp_l$robust_rse <- unname(100 * apply(le, 2, stats::mad)[cmp_l$parameter])
  cmp_l$months <- L
  check[[length(check) + 1L]] <- cmp_l
  cat(sprintf("\n%d months: %d of 1000 fits converged\n", L, cmp_l$converged[1]))
  print(cmp_l[, c("parameter", "expected_rse", "empirical_rse", "robust_rse", "ratio")], row.names = FALSE, digits = 3)
}
utils::write.csv(do.call(rbind, check), file.path(out_dir, "psa-simulation-check.csv"), row.names = FALSE)
