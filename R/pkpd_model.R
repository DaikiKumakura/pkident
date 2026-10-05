#' Define a model for identifiability analysis
#'
#' Reads an ODE model through rxode2 and records, as symbolic expressions, the
#' state equations, the observed outputs, the initial conditions, the known
#' doses and input signals, and the unknown parameters. No identifiability
#' decision is made here.
#'
#' @param model An rxode2 model: a compiled `rxode2` object, model code as a
#'   character string, or a model function with `ini()` and `model()` blocks
#'   (residual error lines containing `~` are ignored).
#' @param outputs Names of the observed quantities: variables assigned in the
#'   model or state names.
#' @param doses A list of [bolus_dose()] and [infusion_dose()] specifications.
#' @param inputs Names of known input signals (for example a plasma
#'   concentration `Cp` driving a pharmacodynamic model). Their values and time
#'   derivatives are treated as known.
#' @param known Names of known constants (for example `DOSE`).
#' @param parameters Names of the unknown parameters to analyse. By default,
#'   every symbol that is not a state, an input, a known constant or time.
#' @param initial A named list of initial conditions as character expressions,
#'   for example `list(E = "kin/kout")`. States without a dose or an initial
#'   condition start at zero.
#' @return An object of class `pkident_model`.
#' @examples
#' m <- pkpd_model(
#'   "d/dt(depot) = -ka*depot
#'    d/dt(central) = ka*depot - CL/V*central
#'    cp = central/V",
#'   outputs = "cp",
#'   doses = list(bolus_dose("depot", "F*DOSE")),
#'   known = "DOSE"
#' )
#' m
#' @export
pkpd_model <- function(model, outputs, doses = list(), inputs = character(),
                       known = character(), parameters = NULL, initial = list()) {
  if (inherits(doses, "pkident_dose")) doses <- list(doses)
  if (!is.character(outputs) || !length(outputs)) {
    pki_abort("PKI009", "`outputs` must name at least one observed quantity.")
  }
  if (!all(vapply(doses, inherits, logical(1), "pkident_dose"))) {
    pki_abort("PKI009", "`doses` must be a list of bolus_dose() or infusion_dose() specifications.")
  }
  if (length(initial) && (is.null(names(initial)) || any(!nzchar(names(initial))))) {
    pki_abort("PKI009", "`initial` must be a named list, for example list(E = \"kin/kout\").")
  }

  code <- model_code(model)
  rx <- tryCatch(
    suppressMessages(rxode2::rxode2(code)),
    error = function(e) pki_abort("PKI001", paste("rxode2 could not parse the model:", conditionMessage(e)))
  )
  mv <- rxode2::rxModelVars(rx)
  states <- as.character(mv$state)
  lhs <- as.character(mv$lhs)
  sym <- suppressMessages(rxode2::rxS(rx))

  # State equations and initial conditions from the model itself
  odes <- stats::setNames(lapply(states, function(s) get(paste0("rx__d_dt_", s, "__"), envir = sym)), states)
  init <- stats::setNames(lapply(states, function(s) {
    nm <- paste0("rx_", s, "_ini_0__")
    if (exists(nm, envir = sym, inherits = FALSE)) get(nm, envir = sym) else symengine::S("0")
  }), states)
  init_from_model <- stats::setNames(
    vapply(states, function(s) exists(paste0("rx_", s, "_ini_0__"), envir = sym, inherits = FALSE), logical(1)),
    states
  )

  # Outputs
  missing_out <- setdiff(outputs, c(lhs, states))
  if (length(missing_out)) {
    pki_abort("PKI002", sprintf("Output(s) not assigned in the model and not states: %s.",
                                paste(missing_out, collapse = ", ")))
  }
  out <- stats::setNames(lapply(outputs, function(o) {
    if (o %in% states) symengine::S(o) else get(o, envir = sym)
  }), outputs)

  # Equations before doses, used when doses are given as event tables
  base_odes <- odes
  base_init <- init

  # Doses
  dose_set <- character()
  for (d in doses) {
    if (!d$state %in% states) {
      pki_abort("PKI004", sprintf("Dose state `%s` is not a state of the model (states: %s).",
                                  d$state, paste(states, collapse = ", ")))
    }
    expr <- parse_expr(d$expression, sprintf("dose into `%s`", d$state))
    if (d$type == "bolus") {
      if (init_from_model[[d$state]] || d$state %in% dose_set) {
        pki_abort("PKI003", sprintf("State `%s` has more than one initial condition (model, bolus or initial).", d$state))
      }
      init[[d$state]] <- expr
      dose_set <- c(dose_set, d$state)
    } else {
      odes[[d$state]] <- odes[[d$state]] + expr
    }
  }

  # Explicit initial conditions
  for (s in names(initial)) {
    if (!s %in% states) {
      pki_abort("PKI004", sprintf("Initial condition given for `%s`, which is not a state.", s))
    }
    if (init_from_model[[s]] || s %in% dose_set) {
      pki_abort("PKI003", sprintf("State `%s` has more than one initial condition (model, bolus or initial).", s))
    }
    init[[s]] <- parse_expr(initial[[s]], sprintf("initial condition of `%s`", s))
    base_init[[s]] <- init[[s]]
    dose_set <- c(dose_set, s)
  }

  # Symbols
  time_symbols <- c("t", "time")
  all_symbols <- unique(c(
    unlist(lapply(c(odes, out, init), symbols_of)),
    character()
  ))
  free <- setdiff(all_symbols, c(states, time_symbols))
  bad_input <- setdiff(inputs, free)
  if (length(bad_input)) {
    pki_abort("PKI005", sprintf("Input signal(s) not used in the model: %s.", paste(bad_input, collapse = ", ")))
  }
  bad_known <- setdiff(known, free)
  if (length(bad_known)) {
    pki_abort("PKI007", sprintf("Known constant(s) not used in the model: %s.", paste(bad_known, collapse = ", ")))
  }
  candidate <- setdiff(free, c(inputs, known))
  if (is.null(parameters)) {
    parameters <- sort(candidate, method = "radix")  # locale-independent order
  } else {
    extra <- setdiff(parameters, candidate)
    if (length(extra)) {
      pki_abort("PKI009", sprintf("Parameter(s) not free symbols of the model: %s.", paste(extra, collapse = ", ")))
    }
    unlisted <- setdiff(candidate, parameters)
    if (length(unlisted)) {
      pki_abort("PKI009", sprintf(paste0("Symbol(s) %s are neither parameters, inputs nor known constants. ",
                                         "Declare them in `known` or `parameters`."), paste(unlisted, collapse = ", ")))
    }
  }
  if (!length(parameters)) pki_abort("PKI006", "The model has no unknown parameters.")

  for (o in names(out)) {
    used <- symbols_of(out[[o]])
    if (o %in% states) next
    reach <- c(states, inputs)
    if (!any(used %in% reach)) {
      pki_abort("PKI008", sprintf("Output `%s` depends on neither a state nor an input signal.", o))
    }
  }

  structure(
    list(
      states = states,
      odes = odes,
      outputs = out,
      initial = init,
      parameters = parameters,
      inputs = inputs,
      known = known,
      doses = doses,
      time_dependent = any(time_symbols %in% all_symbols),
      base = list(odes = base_odes, initial = base_init),
      source = list(code = code, rxode2_version = as.character(utils::packageVersion("rxode2")),
                    symengine_version = as.character(utils::packageVersion("symengine")))
    ),
    class = "pkident_model"
  )
}

#' @export
print.pkident_model <- function(x, ...) {
  cat("<pkident_model>\n")
  cat("States:     ", if (length(x$states)) paste(x$states, collapse = ", ") else "(none)", "\n")
  cat("Outputs:    ", paste(names(x$outputs), collapse = ", "), "\n")
  cat("Parameters: ", paste(x$parameters, collapse = ", "), "\n")
  if (length(x$inputs)) cat("Inputs:     ", paste(x$inputs, collapse = ", "), "(known signals)\n")
  if (length(x$known)) cat("Known:      ", paste(x$known, collapse = ", "), "\n")
  for (s in x$states) cat(sprintf("  d/dt(%s) = %s;  %s(0) = %s\n", s, as.character(x$odes[[s]]), s, as.character(x$initial[[s]])))
  for (o in names(x$outputs)) cat(sprintf("  %s = %s\n", o, as.character(x$outputs[[o]])))
  invisible(x)
}

# Helpers ------------------------------------------------------------------

model_code <- function(model) {
  if (inherits(model, "rxode2")) return(as.character(rxode2::rxNorm(model)))
  if (is.character(model)) return(paste(model, collapse = "\n"))
  if (is.function(model)) {
    ui <- tryCatch(suppressMessages(rxode2::rxode2(model)),
                   error = function(e) pki_abort("PKI001", paste("rxode2 could not read the model function:", conditionMessage(e))))
    lines <- vapply(ui$lstExpr, function(e) paste(deparse(e, width.cutoff = 500L), collapse = " "), character(1))
    return(paste(lines[!grepl("~", lines, fixed = TRUE)], collapse = "\n"))
  }
  pki_abort("PKI009", "`model` must be an rxode2 object, model code or a model function.")
}

parse_expr <- function(x, what) {
  tryCatch(symengine::S(x), error = function(e) pki_abort("PKI001", sprintf("Could not parse %s: `%s`.", what, x)))
}

symbols_of <- function(expr) {
  vapply(as.list(symengine::free_symbols(expr)), as.character, character(1))
}
