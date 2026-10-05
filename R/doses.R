#' Dose specifications
#'
#' Describe how a known dose enters the model.
#'
#' * `bolus()` sets the initial amount of a state (an instantaneous dose at
#'   time zero), for example `bolus("depot", "F*DOSE")`.
#' * `infusion()` adds a known zero-order input rate to a state, for example
#'   `infusion("central", "RATE")`.
#'
#' Symbols used in `amount` or `rate` that are not states become parameters
#' unless they are declared as `known` in [pkpd_model()].
#'
#' @param state Name of the state that receives the dose.
#' @param amount Amount at time zero, as a character expression.
#' @param rate Zero-order input rate, as a character expression.
#' @return An object of class `pkident_dose`.
#' @examples
#' bolus("depot", "F*DOSE")
#' infusion("central", "RATE")
#' @export
bolus <- function(state, amount = "DOSE") {
  check_string(state, "state")
  check_string(amount, "amount")
  structure(list(type = "bolus", state = state, expression = amount), class = "pkident_dose")
}

#' @rdname bolus
#' @export
infusion <- function(state, rate = "RATE") {
  check_string(state, "state")
  check_string(rate, "rate")
  structure(list(type = "infusion", state = state, expression = rate), class = "pkident_dose")
}

#' @export
print.pkident_dose <- function(x, ...) {
  cat(sprintf("<pkident dose> %s into %s: %s\n", x$type, x$state, x$expression))
  invisible(x)
}

check_string <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    pki_abort("PKI009", sprintf("`%s` must be a single non-empty string.", what))
  }
}
