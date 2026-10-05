#' Classify likelihood profiles of an nlmixr2 fit
#'
#' Fixes each selected population parameter at values on either side of its
#' estimate, re-estimates all other parameters with
#' `nlmixr2extra::profileFixed()`, and classifies the change in objective
#' function value (dOFV) following Raue et al. (2009, Bioinformatics
#' 25:1923-1929):
#'
#' * `identifiable`: dOFV exceeds `threshold` on both sides. An approximate
#'   profile-likelihood confidence interval is reported (interpolated on the
#'   signed-root scale; use `profile(fit, method = "llp")` from nlmixr2extra
#'   for a precise interval).
#' * `non_identifiable`: dOFV exceeds the threshold on one side only (the
#'   open direction is reported), or on neither side. A profile that stays
#'   within `flat_tol` of its minimum is reported as flat, which suggests
#'   structural non-identifiability; otherwise it is reported as shallow.
#' * `unresolved`: a profile fit failed.
#'
#' Statements are limited to the explored range, which is reported.
#'
#' Every profile fit is an upper bound of the profile likelihood: a failed or
#' stuck re-estimation can only overstate dOFV. Each point more than
#' `flat_tol` above the original fit, or above the lowest OFV found, is
#' therefore re-estimated from the estimates of a neighbouring profile point
#' (a warm start), and the lower value is kept. If profiling
#' finds a lower OFV than the original fit, dOFV is measured from that lower
#' value and the drop is reported (`ofv_drop`): the original fit was not at
#' its minimum, which is common along non-identifiable directions.
#'
#' By default, each side is explored by stepping away from the estimate by
#' `steps` on the estimation scale, stopping once dOFV exceeds the threshold.
#' The default steps suit log-transformed parameters (`cl <- exp(tcl + ...)`);
#' for other parameters give `grid`.
#'
#' If a [structural_identifiability()] result is given with `map`, each
#' profile is placed next to the structural decision, which distinguishes
#' "structurally identifiable but not determined by these data" from
#' structural non-identifiability.
#'
#' @param fit An nlmixr2 fit (FOCEi or another estimation method supported by
#'   `nlmixr2extra::profileFixed()`).
#' @param which Names of population parameters (thetas) to profile. By default
#'   all estimated thetas that are not residual error parameters.
#' @param grid Optional named list of values (estimation scale) per
#'   parameter; replaces the stepwise search for that parameter.
#' @param threshold dOFV threshold; 3.84 is the 95% point of a chi-squared
#'   distribution with one degree of freedom.
#' @param steps Distances from the estimate, explored on each side in order.
#' @param flat_tol dOFV below which a profile counts as flat (allows for the
#'   convergence noise of repeated estimation).
#' @param structural Optional result of [structural_identifiability()].
#' @param map Named character vector from fit parameter names to pkident
#'   parameter names, for example `c(tcl = "CL", tv = "V")`.
#' @return An object of class `pkident_profiles`. Use [as_table()] with
#'   `"parameters"` or `"profiles"`.
#' @examples
#' \dontrun{
#' fit <- nlmixr2(one_compartment, data, est = "focei")
#' classify_profiles(fit, which = c("tcl", "tv"))
#' }
#' @export
classify_profiles <- function(fit, which = NULL, grid = NULL, threshold = 3.84,
                              steps = c(0.1, 0.25, 0.5, 1, 2), flat_tol = 1,
                              structural = NULL, map = NULL) {
  if (!requireNamespace("nlmixr2extra", quietly = TRUE) || !requireNamespace("nlmixr2est", quietly = TRUE)) {
    pki_abort("PKI011", "classify_profiles() needs the nlmixr2est and nlmixr2extra packages.")
  }
  ini <- tryCatch(fit$iniDf, error = function(e) NULL)
  ofv0 <- tryCatch(fit$objective, error = function(e) NULL)
  if (is.null(ini) || !is.numeric(ofv0) || !is.finite(ofv0)) {
    pki_abort("PKI009", "`fit` must be an nlmixr2 fit with an objective function value.")
  }
  thetas <- ini[!is.na(ini$ntheta) & !ini$fix & is.na(ini$err), , drop = FALSE]
  if (is.null(which)) which <- thetas$name
  bad <- setdiff(which, thetas$name)
  if (length(bad)) {
    pki_abort("PKI009", sprintf("Not estimated population parameters of the fit: %s (available: %s).",
                                paste(bad, collapse = ", "), paste(thetas$name, collapse = ", ")))
  }
  if (!is.null(grid) && (!is.list(grid) || length(setdiff(names(grid), which)))) {
    pki_abort("PKI009", "`grid` must be a list named by parameters in `which`.")
  }
  if (!is.null(structural) && !inherits(structural, "pkident_structural")) {
    pki_abort("PKI009", "`structural` must be a result of structural_identifiability().")
  }
  if (!is.null(structural) && (is.null(map) || is.null(names(map)))) {
    pki_abort("PKI009", "Give `map` (fit parameter names to pkident parameter names) with `structural`.")
  }
  started <- Sys.time()

  # Cold evaluation: nlmixr2extra::profileFixed() starts from the original
  # estimates. Warm evaluation: the same re-estimation started from the
  # estimates of a neighbouring profile point.
  eval_cold <- function(p, value) {
    w <- stats::setNames(data.frame(value), p)
    out <- NULL
    tryCatch(
      suppressWarnings(suppressMessages(utils::capture.output(
        out <- nlmixr2extra::profileFixed(fit, w), type = "output"
      ))),
      error = function(e) NULL
    )
    if (is.null(out) || !is.finite(out$OFV[1])) return(list(dofv = NA_real_, start = NULL))
    keep <- intersect(setdiff(thetas$name, p), names(out))
    list(dofv = out$OFV[1] - ofv0, start = unlist(out[1, keep, drop = TRUE]))
  }
  eval_warm <- function(p, value, start) {
    if (is.null(start) || !length(start)) return(list(dofv = NA_real_, start = NULL))
    args <- c(list(x = fit), stats::setNames(list(str2lang(sprintf("fixed(%.17g)", value))), p), as.list(start))
    ctl <- fit$control
    ctl$print <- 0L
    ctl$covMethod <- 0L
    ctl$calcTables <- FALSE
    ctl$compress <- FALSE
    nf <- NULL
    tryCatch(
      suppressWarnings(suppressMessages(utils::capture.output({
        ui <- do.call(rxode2::ini, args)
        nf <- nlmixr2est::nlmixr2(ui, est = fit$est, control = ctl)
      }, type = "output"))),
      error = function(e) NULL
    )
    if (is.null(nf) || !is.numeric(nf$objective) || !is.finite(nf$objective)) return(list(dofv = NA_real_, start = NULL))
    fe <- nlmixr2est::fixef(nf)
    list(dofv = nf$objective - ofv0, start = fe[intersect(setdiff(thetas$name, p), names(fe))])
  }

  profiles <- list()
  rows <- list()
  for (p in which) {
    est <- thetas$est[thetas$name == p]
    base_start <- stats::setNames(thetas$est[thetas$name != p], thetas$name[thetas$name != p])
    pts <- data.frame(value = est, dofv = 0, warm_start = FALSE, checked = FALSE)
    starts <- list(base_start)
    add_point <- function(x, prev_start) {
      r <- eval_cold(p, x)
      warm <- FALSE
      checked <- FALSE
      if (is.na(r$dofv) || r$dofv > flat_tol) {
        checked <- TRUE
        w <- eval_warm(p, x, prev_start)
        if (!is.na(w$dofv) && (is.na(r$dofv) || w$dofv < r$dofv)) {
          r <- w
          warm <- TRUE
        }
      }
      pts[nrow(pts) + 1L, ] <<- list(x, r$dofv, warm, checked)
      starts[nrow(pts)] <<- list(r$start)
      r
    }
    if (!is.null(grid[[p]])) {
      prev <- base_start
      for (x in sort(grid[[p]])) {
        r <- add_point(x, prev)
        if (!is.null(r$start)) prev <- r$start
      }
    } else {
      for (side in c(-1, 1)) {
        prev <- base_start
        for (st in steps) {
          r <- add_point(est + side * st, prev)
          if (is.na(r$dofv) || r$dofv > threshold) break
          if (!is.null(r$start)) prev <- r$start
        }
      }
    }
    # Second pass: a point more than flat_tol above the lowest OFV found (a
    # lower OFV may have been found after it was evaluated) is re-estimated
    # from its neighbour with the lower OFV.
    ord <- order(pts$value)
    for (k in seq_along(ord)) {
      i <- ord[k]
      if (pts$checked[i] || is.na(pts$dofv[i]) || pts$dofv[i] - min(pts$dofv, na.rm = TRUE) <= flat_tol) next
      nb <- ord[c(k - 1L, k + 1L)[c(k > 1L, k < length(ord))]]
      nb <- nb[!is.na(pts$dofv[nb]) & !vapply(starts[nb], is.null, logical(1))]
      if (!length(nb)) next
      nb <- nb[which.min(pts$dofv[nb])]
      pts$checked[i] <- TRUE
      w <- eval_warm(p, pts$value[i], starts[[nb]])
      if (!is.na(w$dofv) && w$dofv < pts$dofv[i]) {
        pts$dofv[i] <- w$dofv
        pts$warm_start[i] <- TRUE
        starts[i] <- list(w$start)
      }
    }
    pts <- pts[order(pts$value), , drop = FALSE]
    cls <- classify_profile_points(pts$value, pts$dofv, threshold, flat_tol)
    profiles[[p]] <- data.frame(parameter = p, value = pts$value, dofv = pts$dofv,
                                dofv_from_minimum = pts$dofv + cls$drop, warm_start = pts$warm_start,
                                stringsAsFactors = FALSE)
    rows[[p]] <- data.frame(parameter = p, estimate = est, status = cls$status, shape = cls$shape,
                            open_direction = cls$open, lower = cls$lower, upper = cls$upper,
                            profile_minimum = cls$centre, ofv_drop = cls$drop,
                            explored_from = min(pts$value), explored_to = max(pts$value),
                            n_fits = nrow(pts) - 1L, stringsAsFactors = FALSE)
  }
  params <- do.call(rbind, rows)
  rownames(params) <- NULL

  if (!is.null(structural)) {
    sp <- structural$parameters
    params$structural_parameter <- unname(map[params$parameter])
    params$structural_status <- sp$status[match(params$structural_parameter, sp$parameter)]
    params$interpretation <- unname(mapply(interpret_profile, params$status, params$structural_status))
  }

  structure(
    list(
      parameters = params,
      profiles = do.call(rbind, c(unname(profiles), list(make.row.names = FALSE))),
      objective = ofv0,
      settings = list(threshold = threshold, steps = steps, flat_tol = flat_tol),
      versions = list(pkident = as.character(utils::packageVersion("pkident")),
                      nlmixr2extra = as.character(utils::packageVersion("nlmixr2extra"))),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ),
    class = "pkident_profiles"
  )
}

# Classification of one profile. `value` is on the estimation scale and
# `dofv` is relative to the original fit (0 at the estimate). Each profile fit
# gives an upper bound of the profile, so a lower OFV found while profiling is
# real: dOFV is re-referenced to the lowest value found (`drop` > 0 means the
# original fit was not at its minimum). Interval ends are interpolated on the
# signed-root scale, which is exact for a quadratic profile.
classify_profile_points <- function(value, dofv, threshold, flat_tol) {
  out <- list(status = "unresolved", shape = NA_character_, open = NA_character_,
              lower = NA_real_, upper = NA_real_, centre = NA_real_, drop = 0)
  if (any(is.na(dofv))) {
    out$shape <- "failed fit"
    return(out)
  }
  out$drop <- max(0, -min(dofv))
  d <- dofv + out$drop
  k0 <- which.min(d)
  out$centre <- value[k0]
  n <- length(value)
  crosses_lo <- k0 > 1L && any(d[seq_len(k0 - 1L)] > threshold)
  crosses_hi <- k0 < n && any(d[(k0 + 1L):n] > threshold)
  crossing <- function(idx) {
    r <- sqrt(d)
    for (k in seq_along(idx)[-1]) {
      a <- idx[k - 1L]
      b <- idx[k]
      if (d[b] > threshold && d[a] <= threshold) {
        return(value[a] + (sqrt(threshold) - r[a]) / (r[b] - r[a]) * (value[b] - value[a]))
      }
    }
    NA_real_
  }
  if (crosses_lo) out$lower <- crossing(rev(seq_len(k0)))
  if (crosses_hi) out$upper <- crossing(k0:n)
  if (crosses_lo && crosses_hi) {
    out$status <- "identifiable"
    out$shape <- "bounded"
  } else if (crosses_lo || crosses_hi) {
    out$status <- "non_identifiable"
    out$shape <- "one-sided"
    out$open <- if (crosses_lo) "upper" else "lower"
  } else if (max(d) < flat_tol) {
    out$status <- "non_identifiable"
    out$shape <- "flat"
    out$open <- "both"
  } else {
    out$status <- "non_identifiable"
    out$shape <- "shallow"
    out$open <- "both"
  }
  out
}

interpret_profile <- function(profile_status, structural_status) {
  if (is.na(structural_status)) return(NA_character_)
  if (profile_status == "unresolved" || structural_status == "unresolved") return("unresolved")
  s_ok <- structural_status == "identifiable"
  p_ok <- profile_status == "identifiable"
  if (s_ok && p_ok) return("determined by these data")
  if (s_ok && !p_ok) return("structurally identifiable but not determined by these data")
  if (!s_ok && !p_ok) return("structurally non-identifiable")
  "conflict: profile bounded although structurally non-identifiable (check the mapping and the fit)"
}

#' @export
print.pkident_profiles <- function(x, ...) {
  cat("<pkident likelihood profiles>\n")
  cat(sprintf("Threshold dOFV %.2f (Raue et al. 2009); statements apply to the explored range only.\n",
              x$settings$threshold))
  p <- x$parameters
  for (i in seq_len(nrow(p))) {
    extra <- switch(p$shape[i],
      bounded = sprintf("approx. interval %.3g to %.3g", p$lower[i], p$upper[i]),
      `one-sided` = sprintf("open towards %s values", p$open_direction[i]),
      flat = "flat profile (structural non-identifiability suspected)",
      shallow = "threshold not reached on either side",
      p$shape[i])
    cat(sprintf("  %-10s %-17s %s (explored %.3g to %.3g, %d fits)\n", p$parameter[i], p$status[i], extra,
                p$explored_from[i], p$explored_to[i], p$n_fits[i]))
    if (is.finite(p$ofv_drop[i]) && p$ofv_drop[i] > x$settings$flat_tol) {
      cat(sprintf("  %-10s    profiling found an OFV %.2f lower than the original fit (at %.3g)\n", "",
                  p$ofv_drop[i], p$profile_minimum[i]))
    }
    if (!is.null(p$interpretation) && !is.na(p$interpretation[i])) {
      cat(sprintf("  %-10s -> %s (structural: %s)\n", "", p$interpretation[i], p$structural_status[i]))
    }
  }
  invisible(x)
}

#' @rdname as_table
#' @export
as_table.pkident_profiles <- function(x, what = c("parameters", "profiles"), ...) {
  what <- match.arg(what)
  tibble::as_tibble(x[[what]])
}
