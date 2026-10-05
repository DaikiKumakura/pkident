#' One experiment of a design
#'
#' A group of individuals with the same dosing and sampling.
#'
#' @param times Sampling times: a numeric vector used for every output, or a
#'   named list with one vector per observed output. Outputs not named are
#'   not observed in this experiment.
#' @param events Optional dosing event table (see [practical_identifiability()]).
#' @param values Optional named values that override the design values for
#'   this experiment, typically known constants such as a dose `DOSE`.
#' @param n Number of individuals.
#' @param inputs Optional known input signals for this experiment (see
#'   [practical_identifiability()]).
#' @return An object of class `pkident_experiment`.
#' @examples
#' experiment(c(0.5, 1, 2, 4, 8, 24), values = c(DOSE = 100), n = 12)
#' @export
experiment <- function(times, events = NULL, values = NULL, n = 1, inputs = NULL) {
  if (!(is.numeric(times) || (is.list(times) && !is.null(names(times))))) {
    pki_abort("PKI009", "`times` must be a numeric vector or a list named by outputs.")
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) pki_abort("PKI009", "`n` must be a positive number.")
  if (!is.null(values) && (!is.numeric(values) || is.null(names(values)))) {
    pki_abort("PKI009", "`values` must be a named numeric vector.")
  }
  structure(list(times = times, events = events, values = values, n = n, inputs = inputs),
            class = "pkident_experiment")
}

#' @export
print.pkident_experiment <- function(x, ...) {
  tt <- if (is.numeric(x$times)) list(`all outputs` = x$times) else x$times
  cat(sprintf("<pkident experiment> %s individual(s)%s\n", format(x$n),
              if (is.null(x$events)) "" else sprintf(", %d dosing events", nrow(x$events))))
  for (o in names(tt)) cat(sprintf("  %s: %d samples, %g to %g\n", o, length(tt[[o]]), min(tt[[o]]), max(tt[[o]])))
  invisible(x)
}

#' Compare candidate designs
#'
#' For each candidate design supplied by the user (for example the current
#' design, the same design with an added measured output, with added sampling
#' times, or with an added dose level), reports the structural decision for
#' the outputs it observes and the expected relative standard errors from the
#' Fisher information summed over its experiments. No design is optimised;
#' only the supplied candidates are compared.
#'
#' The structural decision depends only on which outputs are observed (it
#' assumes continuous, noise-free observation after the doses of the model
#' definition). Added sampling times and dose levels therefore change only the
#' expected RSEs.
#'
#' @param m A `pkident_model` from [pkpd_model()] that defines every output
#'   that any candidate observes.
#' @param candidates A named list of candidates; each is an [experiment()] or
#'   a list of experiments.
#' @param values Named values of all parameters and known constants (the
#'   design point).
#' @param error A [residual_error()] for every output, or a named list with one
#'   per output.
#' @param inputs Known input signals shared by all experiments.
#' @param structural If `TRUE`, run [structural_identifiability()] for each
#'   distinct set of observed outputs.
#' @param method Method passed to [structural_identifiability()].
#' @param rse_limit RSE (%) above which a parameter is `unresolved`.
#' @param seed Seed for the structural analysis.
#' @param rtol,atol Solver tolerances.
#' @return An object of class `pkident_designs`. Use [as_table()] with
#'   `"summary"` or `"parameters"`.
#' @examples
#' m <- pkpd_model(
#'   "d/dt(central) = -CL/V*central
#'    cp = central/V",
#'   outputs = "cp", doses = bolus_dose("central", "DOSE"), known = "DOSE"
#' )
#' compare_designs(
#'   m,
#'   candidates = list(
#'     early = experiment(c(0.25, 0.5, 1, 2), n = 10),
#'     extended = experiment(c(0.25, 0.5, 1, 2, 8, 24), n = 10)
#'   ),
#'   values = c(CL = 2, V = 20, DOSE = 100),
#'   error = residual_error(add = 0.05, prop = 0.1),
#'   structural = FALSE
#' )
#' @export
compare_designs <- function(m, candidates, values, error = residual_error(prop = 0.1), inputs = list(),
                            structural = TRUE, method = c("both", "lie", "sensitivity"), rse_limit = 30,
                            seed = 1L, rtol = 1e-10, atol = 1e-12) {
  if (!inherits(m, "pkident_model")) pki_abort("PKI009", "`m` must be a pkident_model from pkpd_model().")
  method <- match.arg(method)
  started <- Sys.time()
  if (!is.list(candidates) || !length(candidates) || is.null(names(candidates)) || any(!nzchar(names(candidates)))) {
    pki_abort("PKI009", "`candidates` must be a named list of experiment() objects or lists of them.")
  }
  candidates <- lapply(candidates, function(cand) if (inherits(cand, "pkident_experiment")) list(cand) else cand)
  ok <- vapply(candidates, function(cand) is.list(cand) && length(cand) &&
                 all(vapply(cand, inherits, logical(1), "pkident_experiment")), logical(1))
  if (!all(ok)) pki_abort("PKI009", sprintf("Candidate(s) %s are not experiment() objects or lists of them.",
                                            paste(names(candidates)[!ok], collapse = ", ")))
  if (!is.numeric(values) || is.null(names(values))) pki_abort("PKI009", "`values` must be a named numeric vector.")
  theta <- m$parameters
  outs <- names(m$outputs)

  structural_cache <- list()
  summary_rows <- list()
  param_rows <- list()
  for (cn in names(candidates)) {
    fim <- matrix(0, length(theta), length(theta), dimnames = list(theta, theta))
    observed <- character()
    n_obs <- 0
    for (e in candidates[[cn]]) {
      v <- values
      if (!is.null(e$values)) v[names(e$values)] <- e$values
      v <- check_values(m, v)
      inp <- if (is.null(e$inputs)) inputs else e$inputs
      fi <- fisher_information(m, e$times, v, error, inp, e$events, rtol, atol)
      fim <- fim + e$n * fi$fim
      n_obs <- n_obs + e$n * fi$n_observations
      observed <- union(observed, names(fi$design))
    }
    observed <- outs[outs %in% observed]
    sm <- summarise_fim(fim, check_values(m, values)[theta], rse_limit)

    s_status <- rep(NA_character_, length(theta))
    s_decision <- NA_character_
    s_rank <- NA_integer_
    if (structural) {
      key <- paste(observed, collapse = "+")
      if (is.null(structural_cache[[key]])) {
        ms <- m
        ms$outputs <- m$outputs[observed]
        structural_cache[[key]] <- structural_identifiability(ms, method = method, seed = seed)
      }
      r <- structural_cache[[key]]
      s_status <- r$parameters$status[match(theta, r$parameters$parameter)]
      s_decision <- r$decision
      if (r$decision == "decided") s_rank <- length(theta) - ncol(r$null_bases[[1]])
    }
    p <- sm$parameters
    summary_rows[[cn]] <- data.frame(
      candidate = cn, outputs = paste(observed, collapse = ", "), experiments = length(candidates[[cn]]),
      observations = n_obs, structural_decision = s_decision, structural_rank = s_rank,
      n_parameters = length(theta), n_identifiable = sum(p$status == "identifiable"),
      max_rse_percent = if (all(is.finite(p$rse_percent))) max(p$rse_percent) else Inf,
      condition_number = sm$condition_number, stringsAsFactors = FALSE
    )
    param_rows[[cn]] <- data.frame(candidate = cn, parameter = theta, value = p$value,
                                   structural_status = s_status, rse_percent = p$rse_percent,
                                   practical_status = p$status, stringsAsFactors = FALSE)
  }
  summ <- do.call(rbind, summary_rows)
  params <- do.call(rbind, param_rows)
  rownames(summ) <- NULL
  rownames(params) <- NULL
  first <- params[params$candidate == names(candidates)[1], ]
  params$rse_ratio_to_first <- params$rse_percent / first$rse_percent[match(params$parameter, first$parameter)]

  structure(
    list(summary = summ, parameters = params, structural = structural_cache,
         settings = list(rse_limit = rse_limit, method = method, seed = seed, rtol = rtol, atol = atol),
         versions = list(pkident = as.character(utils::packageVersion("pkident")),
                         rxode2 = m$source$rxode2_version, symengine = m$source$symengine_version),
         elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))),
    class = "pkident_designs"
  )
}

#' @export
print.pkident_designs <- function(x, ...) {
  cat("<pkident design comparison>\n")
  cat("Expected RSE (%) from the Fisher information; structural status for the observed outputs (LOCAL).\n")
  p <- x$parameters
  cands <- unique(p$candidate)
  width <- max(10L, nchar(cands))
  cat(sprintf("  %-12s", "parameter"), sprintf(paste0("%", width, "s"), cands), "\n", sep = " ")
  for (par in unique(p$parameter)) {
    cells <- vapply(cands, function(cn) {
      row <- p[p$candidate == cn & p$parameter == par, ]
      r <- if (is.finite(row$rse_percent)) sprintf("%.1f", row$rse_percent) else "Inf"
      tag <- if (!is.na(row$structural_status) && row$structural_status != "identifiable") "*" else ""
      paste0(r, tag)
    }, character(1))
    cat(sprintf("  %-12s", par), sprintf(paste0("%", width, "s"), cells), "\n", sep = " ")
  }
  s <- x$summary
  cat(sprintf("  %-12s", "identifiable"), sprintf(paste0("%", width, "s"),
      sprintf("%d/%d", s$n_identifiable, s$n_parameters)), "\n", sep = " ")
  if (any(!is.na(p$structural_status) & p$structural_status != "identifiable")) {
    cat("  * not structurally identifiable with the outputs of that candidate\n")
  }
  cat(sprintf("Practically identifiable: RSE <= %g%%.\n", x$settings$rse_limit))
  invisible(x)
}

#' @rdname as_table
#' @export
as_table.pkident_designs <- function(x, what = c("summary", "parameters"), ...) {
  what <- match.arg(what)
  tibble::as_tibble(x[[what]])
}
