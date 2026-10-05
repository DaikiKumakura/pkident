#' Residual error model
#'
#' Measurement error of an output: `sd = sqrt(add^2 + (prop * f)^2)`, where `f`
#' is the model prediction (the `combined2` form of nlmixr2). Use `add` alone
#' for additive error and `prop` alone for proportional error.
#'
#' @param add Additive standard deviation (in the units of the output).
#' @param prop Proportional standard deviation (for example 0.1 for 10%).
#' @return An object of class `pkident_residual`.
#' @examples
#' residual_error(prop = 0.1)
#' residual_error(add = 0.05, prop = 0.1)
#' @export
residual_error <- function(add = 0, prop = 0) {
  ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  if (!ok(add) || !ok(prop) || add + prop == 0) {
    pki_abort("PKI009", "`add` and `prop` must be non-negative numbers, and at least one must be positive.")
  }
  structure(list(add = add, prop = prop), class = "pkident_residual")
}

#' @export
print.pkident_residual <- function(x, ...) {
  cat(sprintf("<pkident residual error> add = %g, prop = %g\n", x$add, x$prop))
  invisible(x)
}

#' Practical identifiability for a given design
#'
#' Computes the expected Fisher information matrix of the unknown parameters
#' for one typical individual (or `n` individuals with the same design), given
#' the sampling times, the parameter values and the residual error model. The
#' output sensitivities are obtained by solving the forward sensitivity
#' equations with rxode2. Because the residual standard deviation may depend
#' on the prediction (proportional error), the information includes the term
#' from the variance as well as the term from the mean.
#'
#' The result reports, for each parameter, the expected standard error and
#' relative standard error (RSE, %), the correlation matrix, and the
#' eigen-decomposition of the information matrix on the log-parameter scale:
#' the weakest directions show which combinations of parameters the design
#' determines poorly, and the condition number summarises how poorly.
#'
#' Statuses:
#' * `identifiable`: the RSE is finite and at most `rse_limit`.
#' * `non_identifiable`: the information matrix is singular in a direction
#'   involving the parameter; the design cannot determine it (structurally or
#'   because of too few samples).
#' * `unresolved`: the RSE is finite but above `rse_limit`. The expected
#'   information is a local, linear approximation that is unreliable when
#'   uncertainty is this large; examine likelihood profiles instead.
#'
#' The default limit of 30% comes from the validation by repeated estimation
#' on simulated data (`validation/fim-simulation.R`): with every RSE below
#' about 15%, expected and observed RSEs agreed within 10%; with RSEs of
#' 30-40%, estimates had skewed tails and the expected RSEs underestimated the
#' spread, including that of a parameter whose own expected RSE was small.
#' When any parameter is `unresolved`, the printed result therefore warns that
#' the RSEs of all parameters are approximate.
#'
#' Doses are those in the model (bolus doses at time zero and constant-rate
#' infusions); their amounts and rates are known constants supplied in
#' `values`.
#'
#' @param m A `pkident_model` from [pkpd_model()].
#' @param design Sampling times: a numeric vector used for every output, or a
#'   named list with one vector per output.
#' @param values Named numeric vector with a value for every parameter and
#'   every known constant of the model.
#' @param error A [residual_error()] used for every output, or a named list
#'   with one per observed output.
#' @param inputs For models with known input signals, a named list of
#'   character expressions in time `t`, for example
#'   `list(Cp = "10*exp(-0.2*t)")`.
#' @param n Number of individuals with the same design (the information is
#'   multiplied by `n`).
#' @param rse_limit RSE (%) above which a parameter is `unresolved`.
#' @param rtol,atol Solver tolerances.
#' @return An object of class `pkident_practical`. Use [as_table()] with
#'   `"parameters"`, `"directions"` or `"correlation"`.
#' @examples
#' m <- reference_model("A3")
#' p <- practical_identifiability(
#'   m,
#'   design = c(0.5, 1, 2, 4, 8, 12, 24),
#'   values = c(ka = 1, CL = 2, V = 20, DOSE = 100),
#'   error = residual_error(add = 0.05, prop = 0.1)
#' )
#' p
#' @export
practical_identifiability <- function(m, design, values, error = residual_error(prop = 0.1), inputs = list(),
                                      n = 1, rse_limit = 30, rtol = 1e-10, atol = 1e-12) {
  if (!inherits(m, "pkident_model")) pki_abort("PKI009", "`m` must be a pkident_model from pkpd_model().")
  started <- Sys.time()
  theta <- m$parameters
  outs <- names(m$outputs)

  # Arguments
  if (!is.numeric(values) || is.null(names(values))) pki_abort("PKI009", "`values` must be a named numeric vector.")
  missing_v <- setdiff(c(theta, m$known), names(values))
  if (length(missing_v)) pki_abort("PKI009", sprintf("`values` lacks: %s.", paste(missing_v, collapse = ", ")))
  values <- values[c(theta, m$known)]
  if (any(!is.finite(values)) || any(values[theta] == 0)) {
    pki_abort("PKI009", "Parameter values must be finite and non-zero.")
  }
  if (is.numeric(design)) design <- stats::setNames(rep(list(design), length(outs)), outs)
  if (!is.list(design) || is.null(names(design)) || length(setdiff(names(design), outs))) {
    pki_abort("PKI009", sprintf("`design` must be a numeric vector or a list named by outputs (%s).", paste(outs, collapse = ", ")))
  }
  for (o in names(design)) {
    tt <- design[[o]]
    if (!is.numeric(tt) || !length(tt) || any(!is.finite(tt)) || any(tt < 0)) {
      pki_abort("PKI009", sprintf("Sampling times for `%s` must be finite and non-negative.", o))
    }
  }
  if (inherits(error, "pkident_residual")) error <- stats::setNames(rep(list(error), length(design)), names(design))
  if (!is.list(error) || length(setdiff(names(design), names(error))) ||
      !all(vapply(error[names(design)], inherits, logical(1), "pkident_residual"))) {
    pki_abort("PKI009", "`error` must be a residual_error() or a list of them named by the sampled outputs.")
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) pki_abort("PKI009", "`n` must be a positive number.")
  missing_in <- setdiff(m$inputs, names(inputs))
  if (length(missing_in)) {
    pki_abort("PKI009", sprintf("Give the known input signal(s) %s in `inputs` as expressions in t.", paste(missing_in, collapse = ", ")))
  }

  # Sensitivities at the union of sampling times
  times <- sort(unique(unlist(design)))
  input_lines <- vapply(m$inputs, function(u) paste0(u, " = ", inputs[[u]]), character(1))
  input_values <- stats::setNames(lapply(m$inputs, function(u) {
    val <- tryCatch(eval(parse(text = inputs[[u]]), list(t = times), baseenv()),
                    error = function(e) pki_abort("PKI009", sprintf("Input `%s` could not be evaluated: %s", u, conditionMessage(e))))
    rep_len(as.numeric(val), length(times))
  }), m$inputs)
  model <- compile_sensitivity_model(m, input_lines)
  osf <- output_sensitivity_function(m)
  sol <- solve_output_sensitivities(m, model, osf, values, times, input_values, rtol, atol)

  # Expected Fisher information (natural scale)
  np <- length(theta)
  fim <- matrix(0, np, np, dimnames = list(theta, theta))
  n_obs <- 0L
  for (o in names(design)) {
    idx <- match(design[[o]], times)
    f <- sol[[o]]$y[idx]
    S <- sol[[o]]$S[idx, , drop = FALSE]
    e <- error[[o]]
    v <- e$add^2 + (e$prop * f)^2
    if (any(v <= 0)) {
      pki_abort("PKI009", sprintf("The residual variance of `%s` is zero at some sampling times; add an additive component.", o))
    }
    w_mean <- 1 / v                                # from the mean
    w_var <- 2 * e$prop^4 * f^2 / v^2              # from the variance, d(var)/dtheta = 2 prop^2 f S
    fim <- fim + crossprod(S * sqrt(w_mean + w_var))
    n_obs <- n_obs + length(idx)
  }
  fim <- n * fim

  # Log-parameter scale: information for log(theta)
  D <- diag(values[theta], np)
  fim_log <- D %*% fim %*% D
  dimnames(fim_log) <- list(theta, theta)
  eg <- eigen((fim_log + t(fim_log)) / 2, symmetric = TRUE)
  ev <- pmax(eg$values, 0)
  rel <- ev / max(ev[1], .Machine$double.xmin)
  null_dir <- rel < 1e-10
  cond <- if (any(null_dir)) Inf else ev[1] / ev[np]

  # Covariance on the log scale (pseudo-inverse over the non-null directions)
  V <- eg$vectors
  keep <- !null_dir
  cov_log <- V[, keep, drop = FALSE] %*% diag(1 / ev[keep], sum(keep)) %*% t(V[, keep, drop = FALSE])
  dimnames(cov_log) <- list(theta, theta)
  null_comp <- if (any(null_dir)) sqrt(rowSums(V[, null_dir, drop = FALSE]^2)) else rep(0, np)
  in_null <- null_comp > 1e-6
  rse <- 100 * sqrt(pmax(diag(cov_log), 0))
  rse[in_null] <- Inf
  se <- abs(values[theta]) * rse / 100
  status <- ifelse(in_null, "non_identifiable", ifelse(rse <= rse_limit, "identifiable", "unresolved"))
  sdv <- sqrt(diag(cov_log))
  cor <- cov_log / outer(sdv, sdv)
  cor[in_null, ] <- NA_real_
  cor[, in_null] <- NA_real_

  params <- data.frame(parameter = theta, value = unname(values[theta]), se = unname(se),
                       rse_percent = unname(rse), status = unname(status), stringsAsFactors = FALSE)
  directions <- data.frame(
    direction = seq_len(np),
    relative_eigenvalue = rel,
    expected_sd_log = ifelse(null_dir, Inf, 1 / sqrt(ev)),
    loadings = vapply(seq_len(np), function(j) loading_string(V[, j], theta), character(1)),
    stringsAsFactors = FALSE
  )
  directions <- directions[rev(seq_len(np)), , drop = FALSE]   # weakest first
  directions$direction <- seq_len(np)
  rownames(directions) <- NULL

  structure(
    list(
      parameters = params,
      directions = directions,
      correlation = cor,
      fim = fim,
      fim_log = fim_log,
      condition_number = cond,
      n_observations = n_obs,
      design = design,
      values = values,
      error = error,
      settings = list(n = n, rse_limit = rse_limit, rtol = rtol, atol = atol),
      model = m,
      versions = list(pkident = as.character(utils::packageVersion("pkident")),
                      rxode2 = m$source$rxode2_version, symengine = m$source$symengine_version),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ),
    class = "pkident_practical"
  )
}

#' @export
print.pkident_practical <- function(x, ...) {
  cat("<pkident practical identifiability>\n")
  cat(sprintf("Expected Fisher information for %s individual(s), %d observations each; LOCAL (linear) approximation.\n",
              format(x$settings$n), x$n_observations))
  cat(sprintf("Condition number (log scale): %s\n", if (is.finite(x$condition_number)) format(signif(x$condition_number, 3)) else "Inf (singular)"))
  p <- x$parameters
  for (i in seq_len(nrow(p))) {
    cat(sprintf("  %-12s RSE %8s%%  %s\n", p$parameter[i],
                if (is.finite(p$rse_percent[i])) formatC(p$rse_percent[i], format = "f", digits = 1) else "Inf",
                p$status[i]))
  }
  if (any(p$status == "unresolved")) {
    cat("Note: some RSEs exceed ", x$settings$rse_limit, "%; the linear approximation is unreliable and the RSEs of ",
        "all parameters may be underestimated. Examine likelihood profiles.\n", sep = "")
  }
  d <- x$directions[1, ]
  cat(sprintf("Weakest direction: %s (expected SD on log scale %s)\n", d$loadings,
              if (is.finite(d$expected_sd_log)) format(signif(d$expected_sd_log, 3)) else "Inf"))
  invisible(x)
}

#' @rdname as_table
#' @export
as_table.pkident_practical <- function(x, what = c("parameters", "directions", "correlation"), ...) {
  what <- match.arg(what)
  if (what == "correlation") {
    cm <- x$correlation
    return(tibble::as_tibble(data.frame(parameter = rownames(cm), cm, check.names = FALSE, row.names = NULL)))
  }
  tibble::as_tibble(x[[what]])
}

loading_string <- function(v, theta) {
  v <- v / sign(v[which.max(abs(v))])
  big <- order(-abs(v))
  big <- big[abs(v[big]) >= 0.1]
  paste(vapply(big, function(j) sprintf("%+.2f log(%s)", v[j], theta[j]), character(1)), collapse = " ")
}
