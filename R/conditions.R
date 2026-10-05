# Coded conditions -----------------------------------------------------------
#
# Every error raised by pkident carries a code so that scripts and continuous
# integration can react to it without parsing the message.

pki_codes <- c(
  PKI001 = "MODEL_PARSE_FAILED",
  PKI002 = "OUTPUT_NOT_FOUND",
  PKI003 = "INITIAL_CONDITION_CONFLICT",
  PKI004 = "UNKNOWN_STATE",
  PKI005 = "INPUT_NOT_IN_MODEL",
  PKI006 = "NO_PARAMETERS",
  PKI007 = "KNOWN_NOT_IN_MODEL",
  PKI008 = "OUTPUT_WITHOUT_STATE_OR_INPUT",
  PKI009 = "INVALID_ARGUMENT",
  PKI010 = "SENSITIVITY_SOLVE_FAILED",
  PKI011 = "SUGGESTED_PACKAGE_MISSING"
)

pki_abort <- function(code, message, call = NULL) {
  stopifnot(code %in% names(pki_codes))
  cond <- structure(
    class = c(paste0("pkident_", code), "pkident_error", "error", "condition"),
    list(message = paste0("[", code, " ", pki_codes[[code]], "] ", message),
         call = call, code = code, name = pki_codes[[code]])
  )
  stop(cond)
}
