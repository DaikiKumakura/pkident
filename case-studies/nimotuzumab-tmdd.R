# Case study 2 (milestone S7): nimotuzumab target-mediated disposition.
#
# Which parameters of a full target-mediated drug disposition (TMDD) model are
# determined by the public nimoData observations (nlmixr2data; source:
# Rodriguez-Vera et al. 2015), and which added measurement would determine
# the rest?
#
# 1. Fit a full TMDD model to nimoData with nlmixr2 (FOCEi). DV is the log of
#    the free drug concentration; the concentration unit is not documented,
#    so the volume is in dose units per concentration unit.
# 2. Structural identifiability for each set of observed outputs.
# 3. Likelihood profiles of the fit (which parameters these data determine).
# 4. Expected RSEs at the fitted values for the actual nimoData design and
#    for candidate designs with added measurements.
#
# Run from the package root: Rscript case-studies/nimotuzumab-tmdd.R
# Set PKIDENT_SKIP_PROFILES=1 to skip step 3 (about 30-60 minutes).

suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages({library(nlmixr2est); library(rxode2)})
out_dir <- "inst/extdata"
dir.create(out_dir, showWarnings = FALSE)

# 1. Fit ----------------------------------------------------------------------
d <- nlmixr2data::nimoData
d <- d[, c("ID", "TIME", "AMT", "RATE", "DV", "EVID", "DOS")]
d$DV[d$EVID == 1] <- NA

tmdd <- function() {
  ini({
    lcl <- log(0.005); lv <- log(1)
    lkon <- log(0.05); lkoff <- log(0.05); lkdeg <- log(0.05); lr0 <- log(5); lkint <- log(0.05)
    eta.cl ~ 0.1; eta.v ~ 0.1
    add.sd <- 0.2
  })
  model({
    cl <- exp(lcl + eta.cl); v <- exp(lv + eta.v)
    kon <- exp(lkon); koff <- exp(lkoff); kdeg <- exp(lkdeg); r0 <- exp(lr0); kint <- exp(lkint)
    d/dt(A) <- -cl / v * A - v * (kon * A / v * R - koff * P)
    d/dt(R) <- r0 * kdeg - kdeg * R - kon * A / v * R + koff * P
    d/dt(P) <- kon * A / v * R - koff * P - kint * P
    R(0) <- r0
    lc <- log(A / v)
    lc ~ add(add.sd)
  })
}
fit <- suppressWarnings(suppressMessages(nlmixr2(tmdd, d, est = "focei",
                                                 control = foceiControl(print = 0, maxOuterIterations = 2000))))
cat(sprintf("Fit: OFV %.2f\n", fit$objective))
print(fit$parFixedDf[, c("Estimate", "SE", "Back-transformed")])
est <- exp(nlmixr2est::fixef(fit)[c("lcl", "lv", "lkon", "lkoff", "lkdeg", "lr0", "lkint")])
names(est) <- c("cl", "v", "kon", "koff", "kdeg", "r0", "kint")
sd_log <- unname(nlmixr2est::fixef(fit)["add.sd"])
utils::write.csv(data.frame(parameter = rownames(fit$parFixedDf), fit$parFixedDf, check.names = FALSE),
                 file.path(out_dir, "nimotuzumab-fit.csv"), row.names = FALSE)

# 2. Structural identifiability ---------------------------------------------
m <- pkpd_model(
  "d/dt(A) = -cl/v*A - v*(kon*A/v*R - koff*P)
   d/dt(R) = r0*kdeg - kdeg*R - kon*A/v*R + koff*P
   d/dt(P) = kon*A/v*R - koff*P - kint*P
   cp = A/v
   rtot = R + P",
  outputs = c("cp", "rtot", "R", "P"),
  doses = bolus_dose("A", "DOSE"), initial = list(R = "r0"), known = "DOSE"
)
ms <- m
ms$outputs <- m$outputs["cp"]
s_free <- structural_identifiability(ms)
print(s_free)

# 3. Likelihood profiles of the fit -----------------------------------------
map <- c(lcl = "cl", lv = "v", lkon = "kon", lkoff = "koff", lkdeg = "kdeg", lr0 = "r0", lkint = "kint")
if (!identical(Sys.getenv("PKIDENT_SKIP_PROFILES"), "1")) {
  prof <- classify_profiles(fit, which = names(map), structural = s_free, map = map)
  print(prof)
  utils::write.csv(prof$parameters, file.path(out_dir, "nimotuzumab-profiles.csv"), row.names = FALSE)
  utils::write.csv(prof$profiles, file.path(out_dir, "nimotuzumab-profile-points.csv"), row.names = FALSE)
}

# 4. Candidate designs --------------------------------------------------------
# The actual design: each patient's own infusions and sampling times.
subject_experiment <- function(id, outputs = "cp", events_scale = 1) {
  di <- d[d$ID == id, ]
  ev <- data.frame(time = di$TIME[di$EVID == 1], amt = di$AMT[di$EVID == 1] * events_scale,
                   rate = di$RATE[di$EVID == 1] * events_scale, cmt = "A")
  tt <- sort(unique(di$TIME[di$EVID == 0]))
  experiment(stats::setNames(rep(list(tt), length(outputs)), outputs), events = ev)
}
ids <- sort(unique(d$ID))
low_ids <- ids[vapply(ids, function(i) d$DOS[d$ID == i][1] == 50, logical(1))]
current <- lapply(ids, subject_experiment)
with_output <- function(o) lapply(ids, subject_experiment, outputs = c("cp", o))
# An added low-dose group: the 50 mg patients' schedules and sampling at 10 mg.
low_dose <- lapply(low_ids, subject_experiment, events_scale = 10 / 50)

values <- c(est, DOSE = 50)
target_sd <- 0.2   # assumed assay error (log scale) of target measurements
error <- list(cp = residual_error(exp = sd_log), rtot = residual_error(exp = target_sd),
              R = residual_error(exp = target_sd), P = residual_error(exp = target_sd))
cmp <- compare_designs(
  m,
  candidates = list(
    current = current,
    `+ total target` = with_output("rtot"),
    `+ free target` = with_output("R"),
    `+ complex` = with_output("P"),
    `+ 10 mg group` = c(current, low_dose)
  ),
  values = values, error = error
)
print(cmp)
print(as_table(cmp, "summary"))
utils::write.csv(cmp$summary, file.path(out_dir, "nimotuzumab-designs-summary.csv"), row.names = FALSE)
utils::write.csv(cmp$parameters, file.path(out_dir, "nimotuzumab-designs-parameters.csv"), row.names = FALSE)
