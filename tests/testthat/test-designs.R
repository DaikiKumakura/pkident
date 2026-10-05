m1 <- pkpd_model("d/dt(central) = -CL/V*central
                  cp = central/V", outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE")
tt <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
v1 <- c(CL = 2, V = 20, DOSE = 100)

analytic_fim_two_doses <- function(times, d1, d2, t2, prop) {
  k <- 2 / 20
  dose_term <- function(D, t0) {
    dt <- pmax(times - t0, 0); on <- times >= t0
    y <- on * D / 20 * exp(-k * dt)
    cbind(y = y, CL = -dt / 20 * y, V = y * (-1 / 20 + 2 * dt / 400))
  }
  a <- dose_term(d1, 0) + dose_term(d2, t2)
  y <- a[, "y"]; S <- a[, c("CL", "V")]
  v <- (prop * y)^2
  crossprod(S * sqrt(1 / v + 2 * prop^4 * y^2 / v^2))
}

test_that("doses from an event table give the same information as the model dose", {
  e <- residual_error(add = 0.05, prop = 0.1)
  p_model <- practical_identifiability(m1, tt, v1, e)
  p_events <- practical_identifiability(m1, tt, v1, e, events = data.frame(time = 0, amt = 100, cmt = "central"))
  expect_equal(p_events$fim, p_model$fim, tolerance = 1e-6)
})

test_that("repeated doses from an event table agree with the analytic information", {
  ev <- data.frame(time = c(0, 6), amt = c(100, 50), cmt = "central")
  p <- practical_identifiability(m1, tt, v1, residual_error(prop = 0.1), events = ev)
  expect_equal(unname(p$fim), unname(analytic_fim_two_doses(tt, 100, 50, 6, 0.1)), tolerance = 1e-6)
})

test_that("exponential error has no variance term", {
  y <- 100 / 20 * exp(-2 / 20 * tt)
  S <- cbind(CL = -tt / 20 * y, V = y * (-1 / 20 + 2 * tt / 400))
  p <- practical_identifiability(m1, tt, v1, residual_error(exp = 0.2))
  expect_equal(unname(p$fim), unname(crossprod(S / (0.2 * y))), tolerance = 1e-8)
  expect_error(residual_error(exp = 0.2, prop = 0.1), class = "pkident_PKI009")
  expect_output(print(residual_error(exp = 0.2)), "log scale")
})

test_that("event doses that depend on parameters are refused", {
  m2 <- pkpd_model("d/dt(depot) = -ka*depot
                    d/dt(central) = ka*depot - CL/V*central
                    cp = central/V", outputs = "cp", doses = bolus_dose("depot", "F*DOSE"), known = "DOSE")
  expect_error(practical_identifiability(m2, tt, c(ka = 1, CL = 2, V = 20, F = 0.5, DOSE = 100),
                                         events = data.frame(time = 0, amt = 100, cmt = "depot")),
               class = "pkident_PKI009")
  expect_error(practical_identifiability(m1, tt, v1, events = data.frame(time = 0, amt = 100, cmt = "gut")),
               class = "pkident_PKI004")
})

test_that("candidates are compared by structural status and expected RSE", {
  m <- pkpd_model("d/dt(depot) = -ka*depot
                   d/dt(central) = ka*depot - CL/V*central
                   cp = central/V
                   amt_gut = depot", outputs = c("cp", "amt_gut"), doses = bolus_dose("depot", "F*DOSE"), known = "DOSE")
  vals <- c(ka = 1, CL = 2, V = 20, F = 0.7, DOSE = 100)
  r <- compare_designs(
    m,
    candidates = list(
      plasma = experiment(list(cp = tt), n = 10),
      `plasma and gut` = experiment(list(cp = tt, amt_gut = c(0.5, 1, 2)), n = 10),
      `two doses` = list(experiment(list(cp = tt), n = 5), experiment(list(cp = tt), values = c(DOSE = 400), n = 5))
    ),
    values = vals, error = residual_error(add = 0.05, prop = 0.1), method = "lie"
  )
  s <- as_table(r, "summary")
  expect_identical(s$structural_rank, c(3L, 4L, 3L))
  p <- as_table(r, "parameters")
  st <- function(cn, par) p$structural_status[p$candidate == cn & p$parameter == par]
  expect_identical(st("plasma", "F"), "combination_only")
  expect_identical(st("plasma and gut", "F"), "identifiable")
  expect_identical(p$practical_status[p$candidate == "plasma" & p$parameter == "V"], "non_identifiable")
  expect_true(is.finite(p$rse_percent[p$candidate == "plasma and gut" & p$parameter == "V"]))
  # a second dose level does not make F, V or CL identifiable in a linear model
  expect_identical(p$practical_status[p$candidate == "two doses" & p$parameter == "F"], "non_identifiable")
  expect_equal(p$rse_ratio_to_first[p$candidate == "plasma" & p$parameter == "ka"], 1)
  expect_output(print(r), "not structurally identifiable")
  expect_identical(s$observations, c(80, 110, 80))
})

test_that("invalid candidates give coded errors", {
  expect_error(compare_designs(m1, list(experiment(tt)), v1), class = "pkident_PKI009")
  expect_error(compare_designs(m1, list(a = tt), v1), class = "pkident_PKI009")
  expect_error(experiment("x"), class = "pkident_PKI009")
  expect_output(print(experiment(tt, n = 3)), "3 individual")
})
