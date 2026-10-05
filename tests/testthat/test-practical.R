tt <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
vals_a1 <- c(CL = 2, V = 20, DOSE = 100)

analytic_fim_a1 <- function(add, prop) {
  y <- 100 / 20 * exp(-2 / 20 * tt)
  S <- cbind(CL = -tt / 20 * y, V = y * (-1 / 20 + 2 * tt / 400))
  v <- add^2 + (prop * y)^2
  crossprod(S * sqrt(1 / v + 2 * prop^4 * y^2 / v^2))
}

test_that("the Fisher information agrees with the analytic result for one-compartment IV bolus", {
  for (e in list(c(0.05, 0.1), c(0, 0.15), c(0.1, 0))) {
    p <- practical_identifiability(reference_model("A1"), tt, vals_a1, residual_error(add = e[1], prop = e[2]))
    expect_equal(unname(p$fim), unname(analytic_fim_a1(e[1], e[2])), tolerance = 1e-8)
  }
})

test_that("information scales with the number of individuals", {
  e <- residual_error(add = 0.05, prop = 0.1)
  p1 <- practical_identifiability(reference_model("A1"), tt, vals_a1, e)
  p4 <- practical_identifiability(reference_model("A1"), tt, vals_a1, e, n = 4)
  expect_equal(p4$parameters$rse_percent, p1$parameters$rse_percent / 2, tolerance = 1e-10)
})

test_that("expected RSEs agree with the spread of repeated estimates", {
  e <- residual_error(add = 0.05, prop = 0.1)
  m <- reference_model("A3")
  v <- c(ka = 1, CL = 2, V = 20, DOSE = 100)
  p <- practical_identifiability(m, tt, v, e)
  cmp <- compare_rse(p, simulate_and_fit(m, tt, v, e, reps = 300, seed = 11))
  expect_identical(cmp$converged[1], 300L)
  expect_true(all(abs(cmp$ratio - 1) < 0.15), label = paste(round(cmp$ratio, 3), collapse = ", "))
})

test_that("structurally non-identifiable parameters give a singular information matrix", {
  p <- practical_identifiability(reference_model("A2"), tt, c(ka = 1, CL = 2, V = 20, F = 0.7, DOSE = 100))
  st <- stats::setNames(p$parameters$status, p$parameters$parameter)
  expect_identical(unname(st[c("CL", "F", "V")]), rep("non_identifiable", 3))
  expect_identical(unname(st["ka"]), "identifiable")
  expect_identical(p$condition_number, Inf)
  expect_true(all(c("log(CL)", "log(F)", "log(V)") %in% strsplit(gsub("[+-][0-9.]+ ", "", p$directions$loadings[1]), " ")[[1]]))
})

test_that("a sparse design is unresolved, not identifiable", {
  p <- practical_identifiability(reference_model("A3"), c(0.5, 1, 2), c(ka = 1, CL = 2, V = 20, DOSE = 100),
                                 residual_error(add = 0.05, prop = 0.1))
  expect_true(all(p$parameters$status == "unresolved"))
  expect_true(is.finite(p$condition_number) && p$condition_number > 1e4)
  expect_output(print(p), "Examine likelihood profiles")
})

test_that("known input signals and per-output designs are supported", {
  p <- practical_identifiability(reference_model("A10_Rtot_known"), seq(0, 24, 1),
                                 c(ke0 = 0.5, kon = 0.2, koff = 0.5, ke = 1, Rtot = 10),
                                 residual_error(add = 0.02, prop = 0.01), inputs = list(Cp = "10*exp(-0.2*t)"))
  expect_true(all(p$parameters$status == "identifiable"))

  m2 <- pkpd_model("d/dt(central) = -CL/V*central
                    cp = central/V", outputs = c("cp", "central"), doses = bolus_dose("central", "DOSE"), known = "DOSE")
  one <- practical_identifiability(m2, list(cp = tt), vals_a1, residual_error(add = 0.05, prop = 0.1))
  two <- practical_identifiability(m2, list(cp = tt, central = c(1, 4)), vals_a1,
                                   list(cp = residual_error(add = 0.05, prop = 0.1), central = residual_error(prop = 0.1)))
  expect_identical(two$n_observations, 10L)
  expect_true(all(two$parameters$rse_percent < one$parameters$rse_percent))
})

test_that("tables are returned as tibbles", {
  p <- practical_identifiability(reference_model("A1"), tt, vals_a1)
  expect_s3_class(as_table(p, "parameters"), "tbl_df")
  expect_identical(nrow(as_table(p, "directions")), 2L)
  expect_identical(names(as_table(p, "correlation")), c("parameter", "CL", "V"))
  expect_output(print(p), "LOCAL")
})

test_that("invalid arguments give coded errors", {
  m <- reference_model("A1")
  expect_error(residual_error(), class = "pkident_PKI009")
  expect_error(residual_error(add = -1), class = "pkident_PKI009")
  expect_error(practical_identifiability(m, tt, c(CL = 2, V = 20)), class = "pkident_PKI009")
  expect_error(practical_identifiability(m, list(y = tt), vals_a1), class = "pkident_PKI009")
  expect_error(practical_identifiability(m, c(-1, 1), vals_a1), class = "pkident_PKI009")
  expect_error(practical_identifiability(reference_model("A9"), tt, c(ke0 = 1, Kd = 1, ke = 1, Rtot = 1)),
               class = "pkident_PKI009")
})
