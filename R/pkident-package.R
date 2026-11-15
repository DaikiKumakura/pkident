#' pkident: structural and practical identifiability for PK/PD models
#'
#' pkident tells which parameters of a pharmacokinetic or pharmacodynamic
#' ODE model, written in rxode2, the data can determine.
#'
#' @section Workflow:
#' 1. Describe the model: [pkpd_model()] with the observed outputs, the doses
#'    ([bolus_dose()], [infusion_dose()]) and the known quantities, or a
#'    [reference_model()].
#' 2. Structure: [structural_identifiability()] decides, for continuous
#'    noise-free observation, which parameters and which combinations (such
#'    as `CL/F`) can be determined at all. Two independent methods must
#'    agree.
#' 3. Design: [practical_identifiability()] gives the expected relative
#'    standard error of each parameter for a sampling design ([experiment()],
#'    [residual_error()]); [compare_designs()] compares candidates.
#' 4. Fit: [classify_profiles()] classifies the likelihood profiles of an
#'    nlmixr2 fit and links them to the structural result.
#'
#' Every decision is local and says so; anything the methods cannot decide is
#' reported as `unresolved`. [as_table()] turns any result into a data frame.
#'
#' @keywords internal
"_PACKAGE"
