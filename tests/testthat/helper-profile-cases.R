# Synthetic nlmixr2 fits with known identifiability, used by the profile
# tests and by validation/profile-classification.R.
#
# Case "oral":        one-compartment oral, F known, rich sampling.
#                     All parameters identifiable.
# Case "oral_f":      the same data fitted with bioavailability F estimated.
#                     ka identifiable; CL, V and F enter only as CL/F and V/F
#                     (structurally non-identifiable; flat profiles).
# Case "emax_low":    direct Emax effect with a known concentration far below
#                     EC50. Structurally identifiable, but these data do not
#                     bound Emax and EC50 from above.

profile_case_data <- function(case, nid = 24, seed = 1) {
  set.seed(seed)
  if (case %in% c("oral", "oral_f")) {
    sim <- rxode2::rxode2({
      ka <- 1 * exp(eka); cl <- 2 * exp(ecl); v <- 20
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
    })
    times <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    ev <- rxode2::et(rxode2::et(rxode2::et(amt = 100, cmt = "depot"), times), id = seq_len(nid))
    s <- suppressMessages(rxode2::rxSolve(sim, ev, returnType = "data.frame", addDosing = TRUE,
      params = data.frame(id = seq_len(nid), eka = stats::rnorm(nid, 0, 0.3), ecl = stats::rnorm(nid, 0, 0.3))))
    obs <- s$evid == 0
    d <- data.frame(ID = s$id, TIME = s$time, EVID = ifelse(obs, 0L, 1L), AMT = ifelse(obs, 0, 100),
                    CMT = ifelse(obs, "central", "depot"), DV = ifelse(obs, s$cp * (1 + 0.1 * stats::rnorm(nrow(s))), NA))
    return(d)
  }
  if (case == "emax_low") {
    times <- seq(0, 12, 1)
    d <- expand.grid(TIME = times, ID = seq_len(nid))
    dose <- rep(c(0.125, 0.25, 0.5), length.out = nid)[d$ID]
    e0 <- 10 * exp(stats::rnorm(nid, 0, 0.2))[d$ID]
    d$CP <- dose * exp(-0.15 * d$TIME)                  # 0.02 to 0.5, far below EC50 = 20
    eff <- e0 + 50 * d$CP / (20 + d$CP)
    d$DV <- eff + 0.3 * stats::rnorm(nrow(d))
    d$EVID <- 0L
    return(d[, c("ID", "TIME", "EVID", "CP", "DV")])
  }
  stop("unknown case")
}

profile_case_model <- function(case) {
  if (case == "oral") {
    return(function() {
      ini({
        tka <- log(1); tcl <- log(2); tv <- log(20)
        eta.ka ~ 0.09; eta.cl ~ 0.09
        prop.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ prop(prop.sd)
      })
    })
  }
  if (case == "oral_f") {
    return(function() {
      ini({
        tka <- log(1); tcl <- log(2); tv <- log(20); tf <- log(0.7)
        eta.ka ~ 0.09; eta.cl ~ 0.09
        prop.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        f(depot) <- exp(tf)
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ prop(prop.sd)
      })
    })
  }
  if (case == "emax_low") {
    return(function() {
      ini({
        te0 <- log(10); temax <- log(50); tec50 <- log(20)
        eta.e0 ~ 0.04
        add.sd <- 0.3
      })
      model({
        e0 <- exp(te0 + eta.e0); emax <- exp(temax); ec50 <- exp(tec50)
        eff <- e0 + emax * CP / (ec50 + CP)
        eff ~ add(add.sd)
      })
    })
  }
  stop("unknown case")
}

fit_profile_case <- function(case, ...) {
  d <- profile_case_data(case, ...)
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    profile_case_model(case), d, est = "focei", control = nlmixr2est::foceiControl(print = 0)
  )))
}
