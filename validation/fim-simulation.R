# Validation of practical_identifiability() (milestone S5).
#
# For each scenario, the expected relative standard errors from the Fisher
# information are compared with the spread of maximum-likelihood estimates
# from repeated simulated data sets. Run from the package root:
#   Rscript validation/fim-simulation.R
# Results are written to validation/fim-simulation-results.csv. robust_rse is
# 100 x the median absolute deviation of the log estimates, which is
# insensitive to skewed tails.

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-fim-simulation.R")

reps <- 1000
scenarios <- list(
  list(id = "A1", label = "One-compartment IV bolus",
       design = c(0.25, 0.5, 1, 2, 4, 8, 12, 24),
       values = c(CL = 2, V = 20, DOSE = 100),
       error = residual_error(add = 0.05, prop = 0.1)),
  list(id = "A3", label = "One-compartment oral, F known",
       design = c(0.25, 0.5, 1, 2, 4, 8, 12, 24),
       values = c(ka = 1, CL = 2, V = 20, DOSE = 100),
       error = residual_error(add = 0.05, prop = 0.1)),
  list(id = "A4", label = "Two-compartment IV bolus",
       design = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 12, 24, 48),
       values = c(CL = 2, V1 = 10, Q = 3, V2 = 30, DOSE = 100),
       error = residual_error(add = 0.02, prop = 0.1)),
  list(id = "A10_Rtot_known", label = "Effect compartment, dynamic binding (known Cp), precise data",
       design = seq(0, 24, 1),
       values = c(ke0 = 0.5, kon = 0.2, koff = 0.5, ke = 1, Rtot = 10),
       error = residual_error(add = 0.02, prop = 0.01),
       inputs = list(Cp = "10*exp(-0.2*t)")),
  list(id = "A10_Rtot_known", label = "Same with noisier data (limitation: RSEs of 30-40%)",
       design = seq(0, 24, 1),
       values = c(ke0 = 0.5, kon = 0.2, koff = 0.5, ke = 1, Rtot = 10),
       error = residual_error(add = 0.1, prop = 0.05),
       inputs = list(Cp = "10*exp(-0.2*t)")),
  list(id = "A5", label = "Michaelis-Menten elimination, low concentrations (expected: unresolved)",
       design = c(0.25, 0.5, 1, 2, 4, 8, 12, 24),
       values = c(VMAX = 10, KM = 5, V = 20, DOSE = 100),
       error = residual_error(add = 0.05, prop = 0.1))
)

rows <- list()
for (s in scenarios) {
  inputs <- if (is.null(s$inputs)) list() else s$inputs
  p <- practical_identifiability(reference_model(s$id), s$design, s$values, s$error, inputs = inputs)
  t0 <- Sys.time()
  sim <- simulate_and_fit(reference_model(s$id), s$design, s$values, s$error, inputs = inputs, reps = reps, seed = 2026)
  cmp <- compare_rse(p, sim)
  cmp$status <- p$parameters$status
  le <- log(sim$estimates[sim$converged, , drop = FALSE])
  cmp$robust_rse <- unname(100 * apply(le, 2, stats::mad)[cmp$parameter])
  cmp$scenario <- s$id
  cmp$label <- s$label
  rows[[length(rows) + 1L]] <- cmp
  cat(sprintf("\n%s: %s (%d of %d fits converged, %.0f s)\n", s$id, s$label, cmp$converged[1], reps,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(cmp[, c("parameter", "status", "expected_rse", "empirical_rse", "robust_rse", "ratio")], row.names = FALSE, digits = 3)
}
# Cross-check of the fitting routine: refit the first data sets of A3 with
# Nelder-Mead (no derivatives) and compare the estimates.
s <- scenarios[[2]]
m <- reference_model(s$id)
sim <- simulate_and_fit(m, s$design, s$values, s$error, reps = 5, seed = 7)
model <- compile_sensitivity_model(m, character())
osf <- output_sensitivity_function(m)
set.seed(7)
truth <- solve_output_sensitivities(m, model, osf, s$values, s$design, list(), 1e-10, 1e-12)$cp$y
sdv <- sqrt(s$error$add^2 + (s$error$prop * truth)^2)
obs <- matrix(truth, 5, length(truth), byrow = TRUE) + matrix(stats::rnorm(5 * length(truth)), 5) * matrix(sdv, 5, length(truth), byrow = TRUE)
dev <- vapply(1:5, function(i) {
  nll <- function(lt) {
    f <- solve_output_sensitivities(m, model, osf, c(stats::setNames(exp(lt), m$parameters), DOSE = 100),
                                    s$design, list(), 1e-10, 1e-12)$cp$y
    v <- s$error$add^2 + (s$error$prop * f)^2
    sum(log(v) + (obs[i, ] - f)^2 / v)
  }
  nm <- stats::optim(log(sim$estimates[i, ]) + 0.05, nll, control = list(reltol = 1e-14, maxit = 5000))
  max(abs(exp(nm$par) / sim$estimates[i, ] - 1))
}, numeric(1))
cat(sprintf("\nFitting cross-check (A3, 5 data sets): largest relative difference from Nelder-Mead %.1e\n", max(dev)))

res <- do.call(rbind, rows)
utils::write.csv(res[, c("scenario", "label", "parameter", "status", "expected_rse", "empirical_rse", "robust_rse", "ratio", "converged", "reps")],
                 "validation/fim-simulation-results.csv", row.names = FALSE)
