# Validation of classify_profiles() (milestone S6).
#
# Three synthetic nlmixr2 fits whose identifiability is known by construction
# (see tests/testthat/helper-profile-cases.R) are profiled and classified,
# and the classification is compared with the expected result. The profiles
# are also joined with structural_identifiability() results. Run from the
# package root:
#   Rscript validation/profile-classification.R
# Results are written to validation/profile-classification-results.csv and
# validation/profile-classification-profiles.csv.

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-profile-cases.R")

emax_model <- pkpd_model("effect = E0 + Emax*Cp/(EC50 + Cp)", outputs = "effect", inputs = "Cp")

cases <- list(
  oral = list(
    expected = c(tka = "identifiable", tcl = "identifiable", tv = "identifiable"),
    structural = structural_identifiability(reference_model("A3")),
    map = c(tka = "ka", tcl = "CL", tv = "V"),
    truth = c(tka = log(1), tcl = log(2), tv = log(20))
  ),
  oral_f = list(
    expected = c(tka = "identifiable", tcl = "non_identifiable", tv = "non_identifiable", tf = "non_identifiable"),
    expected_shape = c(tcl = "flat", tv = "flat", tf = "flat"),
    structural = structural_identifiability(reference_model("A2")),
    map = c(tka = "ka", tcl = "CL", tv = "V", tf = "F")
  ),
  emax_low = list(
    expected = c(te0 = "identifiable", temax = "non_identifiable", tec50 = "non_identifiable"),
    expected_open = c(temax = "upper", tec50 = "upper"),
    structural = structural_identifiability(emax_model),
    map = c(te0 = "E0", temax = "Emax", tec50 = "EC50")
  )
)

rows <- list()
profiles <- list()
for (case in names(cases)) {
  cs <- cases[[case]]
  t0 <- Sys.time()
  fit <- fit_profile_case(case)
  res <- classify_profiles(fit, which = names(cs$expected), structural = cs$structural, map = cs$map)
  cat(sprintf("\n== %s (OFV %.2f; %.0f s)\n", case, fit$objective, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(res)
  p <- res$parameters
  p$case <- case
  p$expected <- unname(cs$expected[p$parameter])
  p$pass <- p$status == p$expected
  if (!is.null(cs$expected_shape)) {
    k <- p$parameter %in% names(cs$expected_shape)
    p$pass[k] <- p$pass[k] & p$shape[k] == cs$expected_shape[p$parameter[k]]
  }
  if (!is.null(cs$expected_open)) {
    k <- p$parameter %in% names(cs$expected_open)
    p$pass[k] <- p$pass[k] & p$open_direction[k] == cs$expected_open[p$parameter[k]]
  }
  if (!is.null(cs$truth)) {
    p$truth <- unname(cs$truth[p$parameter])
    p$interval_covers_truth <- p$lower <= p$truth & p$truth <= p$upper
  }
  rows[[case]] <- p
  pr <- res$profiles
  pr$case <- case
  profiles[[case]] <- pr
}
cols <- c("case", "parameter", "estimate", "status", "expected", "pass", "shape", "open_direction", "lower", "upper", "ofv_drop", "profile_minimum",
          "explored_from", "explored_to", "n_fits", "structural_status", "interpretation")
out <- do.call(rbind, lapply(rows, function(p) p[, cols]))
cat("\nAll classifications as expected:", all(out$pass), "\n")
cat("Profile intervals cover the true value (case oral):", paste(rows$oral$interval_covers_truth, collapse = ", "), "\n")
utils::write.csv(out, "validation/profile-classification-results.csv", row.names = FALSE)
utils::write.csv(do.call(rbind, c(unname(profiles), list(make.row.names = FALSE)))[, c("case", "parameter", "value", "dofv", "dofv_from_minimum", "warm_start")],
                 "validation/profile-classification-profiles.csv", row.names = FALSE)
