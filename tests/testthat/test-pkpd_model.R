test_that("all reference models load with the expected structure", {
  ids <- reference_models()$id
  expect_length(ids, 11)
  for (id in ids) {
    m <- reference_model(id)
    expect_s3_class(m, "pkident_model")
    exp <- attr(m, "expected")
    expect_setequal(m$parameters, c(exp$identifiable, exp$non_identifiable))
  }
})

test_that("states, outputs, doses and initial conditions are extracted", {
  m <- reference_model("A2")
  expect_identical(m$states, c("depot", "central"))
  expect_identical(names(m$outputs), "cp")
  expect_identical(as.character(symengine::expand(m$initial$depot - symengine::S("F*DOSE"))), "0")
  expect_identical(as.character(m$initial$central), "0")
  expect_identical(m$known, "DOSE")
  expect_setequal(m$parameters, c("CL", "F", "V", "ka"))
  expect_false(m$time_dependent)
})

test_that("intermediate variables are substituted into the equations", {
  m <- pkpd_model("k = CL/V
                   d/dt(central) = -k*central
                   cp = central/V",
                  outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE")
  expect_false("k" %in% m$parameters)
  expect_setequal(m$parameters, c("CL", "V"))
})

test_that("known input signals are separated from parameters", {
  m <- reference_model("A9")
  expect_identical(m$inputs, "Cp")
  expect_false("Cp" %in% m$parameters)
  expect_setequal(m$parameters, c("Kd", "Rtot", "ke", "ke0"))
})

test_that("a state can be observed directly and initial conditions can be given", {
  m <- reference_model("A8")
  expect_identical(as.character(m$outputs$effect), "effect")
  expect_identical(as.character(m$initial$effect), "kin/kout")
})

test_that("models without states are accepted", {
  m <- reference_model("A7")
  expect_length(m$states, 0)
  expect_setequal(m$parameters, c("Kd", "Rtot", "ke"))
})

test_that("an infusion adds its rate to the state equation", {
  m <- pkpd_model("d/dt(central) = -CL/V*central
                   cp = central/V",
                  outputs = "cp", doses = infusion_dose("central", "RATE"), known = "RATE")
  expect_match(as.character(m$odes$central), "RATE", fixed = TRUE)
  expect_identical(as.character(m$initial$central), "0")
})

test_that("compiled rxode2 objects and model functions are accepted", {
  rx <- suppressMessages(rxode2::rxode2("d/dt(central) = -CL/V*central
                                         cp = central/V"))
  m1 <- pkpd_model(rx, outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE")
  expect_setequal(m1$parameters, c("CL", "V"))

  f <- function() {
    ini({
      tcl <- log(1); tv <- log(10); add.sd <- 0.1
    })
    model({
      CL <- exp(tcl); V <- exp(tv)
      d/dt(central) <- -CL/V*central
      cp <- central/V
      cp ~ add(add.sd)
    })
  }
  m2 <- pkpd_model(f, outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE")
  expect_setequal(m2$parameters, c("tcl", "tv"))
})

test_that("errors carry codes", {
  code <- "d/dt(central) = -CL/V*central
           cp = central/V"
  expect_error(pkpd_model(code, outputs = "conc", doses = bolus_dose("central", "DOSE"), known = "DOSE"),
               class = "pkident_PKI002")
  expect_error(pkpd_model(code, outputs = "cp", doses = bolus_dose("gut", "DOSE"), known = "DOSE"),
               class = "pkident_PKI004")
  expect_error(pkpd_model(code, outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE", inputs = "Cp"),
               class = "pkident_PKI005")
  expect_error(pkpd_model(code, outputs = "cp", doses = bolus_dose("central", "DOSE"), known = c("DOSE", "WT")),
               class = "pkident_PKI007")
  expect_error(pkpd_model(code, outputs = "cp", doses = list(bolus_dose("central", "DOSE"), bolus_dose("central", "DOSE")),
                          known = "DOSE"),
               class = "pkident_PKI003")
  expect_error(pkpd_model(code, outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE",
                          parameters = "CL"),
               class = "pkident_PKI009")
  expect_error(pkpd_model("d/dt(central) = -CL/V*central +", outputs = "cp"), class = "pkident_PKI001")
  err <- tryCatch(pkpd_model(code, outputs = "conc"), error = function(e) e)
  expect_identical(err$code, "PKI002")
})
