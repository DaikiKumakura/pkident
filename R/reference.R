#' Reference models with known identifiability
#'
#' The acceptance models used to validate pkident. Each has an expected result
#' fixed from a published table or an analytical derivation before any
#' identifiability code was written (see the package's validation record).
#'
#' `reference_models()` lists them; `reference_model()` builds one as a
#' [pkpd_model()] with the expected result attached as attribute `expected`.
#'
#' @param id Model identifier, for example `"A2"`.
#' @return `reference_models()` returns a data frame; `reference_model()`
#'   returns a `pkident_model`.
#' @references
#' Janzén DLI, et al. Parameter identifiability of fundamental pharmacodynamic
#' models. Front Physiol. 2016;7:590. \doi{10.3389/fphys.2016.00590}
#'
#' Cheung SYA, Yates JWT, Aarons L. The design and analysis of parallel
#' experiments to produce structurally identifiable models. J Pharmacokinet
#' Pharmacodyn. 2013;40:93-100. \doi{10.1007/s10928-012-9291-z}
#' @examples
#' reference_models()
#' reference_model("A2")
#' @export
reference_models <- function() {
  data.frame(
    id = names(reference_specs),
    title = vapply(reference_specs, `[[`, character(1), "title"),
    source = vapply(reference_specs, `[[`, character(1), "source"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname reference_models
#' @export
reference_model <- function(id) {
  check_string(id, "id")
  if (!id %in% names(reference_specs)) {
    pki_abort("PKI009", sprintf("Unknown reference model `%s`. Available: %s.", id, paste(names(reference_specs), collapse = ", ")))
  }
  s <- reference_specs[[id]]
  m <- pkpd_model(s$code, outputs = s$outputs, doses = s$doses, inputs = s$inputs,
                  known = s$known, initial = s$initial)
  attr(m, "expected") <- s$expected
  attr(m, "reference_id") <- id
  m
}

expected <- function(identifiable = character(), non_identifiable = character(),
                     combinations = character(), global_note = NA_character_) {
  list(identifiable = identifiable, non_identifiable = non_identifiable,
       combinations = combinations, global_note = global_note)
}

reference_specs <- list(
  A1 = list(
    title = "One-compartment IV bolus",
    source = "Analytical derivation",
    code = "d/dt(central) = -CL/V*central
            cp = central/V",
    outputs = "cp", doses = list(bolus("central", "DOSE")), inputs = character(), known = "DOSE", initial = list(),
    expected = expected(identifiable = c("CL", "V"))
  ),
  A2 = list(
    title = "One-compartment oral, bioavailability unknown",
    source = "Janzen et al. 2016 (citing Cheung et al. 2013); analytical derivation",
    code = "d/dt(depot) = -ka*depot
            d/dt(central) = ka*depot - CL/V*central
            cp = central/V",
    outputs = "cp", doses = list(bolus("depot", "F*DOSE")), inputs = character(), known = "DOSE", initial = list(),
    expected = expected(identifiable = "ka", non_identifiable = c("CL", "F", "V"), combinations = c("V/F", "CL/F"),
                        global_note = "ka is locally identifiable; flip-flop gives a second global solution")
  ),
  A3 = list(
    title = "One-compartment oral, bioavailability known",
    source = "Analytical derivation",
    code = "d/dt(depot) = -ka*depot
            d/dt(central) = ka*depot - CL/V*central
            cp = central/V",
    outputs = "cp", doses = list(bolus("depot", "DOSE")), inputs = character(), known = "DOSE", initial = list(),
    expected = expected(identifiable = c("CL", "V", "ka"),
                        global_note = "locally identifiable; not globally (flip-flop: (ka, k, V) and (k, ka, V*k/ka))")
  ),
  A4 = list(
    title = "Two-compartment IV bolus, central compartment observed",
    source = "Analytical derivation",
    code = "d/dt(central) = -(CL/V1 + Q/V1)*central + Q/V2*peripheral
            d/dt(peripheral) = Q/V1*central - Q/V2*peripheral
            cp = central/V1",
    outputs = "cp", doses = list(bolus("central", "DOSE")), inputs = character(), known = "DOSE", initial = list(),
    expected = expected(identifiable = c("CL", "Q", "V1", "V2"))
  ),
  A5 = list(
    title = "One-compartment Michaelis-Menten elimination, IV bolus",
    source = "Analytical derivation",
    code = "d/dt(central) = -VMAX*(central/V)/(KM + central/V)
            cp = central/V",
    outputs = "cp", doses = list(bolus("central", "DOSE")), inputs = character(), known = "DOSE", initial = list(),
    expected = expected(identifiable = c("KM", "V", "VMAX"))
  ),
  A6 = list(
    title = "One-compartment parallel linear and Michaelis-Menten elimination, IV bolus",
    source = "Analytical derivation",
    code = "d/dt(central) = -CL/V*central - VMAX*(central/V)/(KM + central/V)
            cp = central/V",
    outputs = "cp", doses = list(bolus("central", "DOSE")), inputs = character(), known = "DOSE", initial = list(),
    expected = expected(identifiable = c("CL", "KM", "V", "VMAX"))
  ),
  A7 = list(
    title = "Direct effect, steady-state receptor binding, linear transduction (known Cp)",
    source = "Janzen et al. 2016, Model 1",
    code = "effect = ke*Rtot*Cp/(Kd + Cp)",
    outputs = "effect", doses = list(), inputs = "Cp", known = character(), initial = list(),
    expected = expected(identifiable = "Kd", non_identifiable = c("Rtot", "ke"), combinations = "Rtot*ke")
  ),
  A8 = list(
    title = "Indirect stimulation of production, steady-state receptor binding (known Cp)",
    source = "Janzen et al. 2016, Model 3 (baseline kin/kout)",
    code = "d/dt(effect) = kin*(1 + ke*Rtot*Cp/(Kd + Cp)) - kout*effect",
    outputs = "effect", doses = list(), inputs = "Cp", known = character(), initial = list(effect = "kin/kout"),
    expected = expected(identifiable = c("Kd", "kin", "kout"), non_identifiable = c("Rtot", "ke"), combinations = "Rtot*ke")
  ),
  A9 = list(
    title = "Effect compartment, steady-state receptor binding, linear transduction (known Cp)",
    source = "Janzen et al. 2016, Model 9",
    code = "d/dt(Ce) = ke0*(Cp - Ce)
            effect = ke*Rtot*Ce/(Kd + Ce)",
    outputs = "effect", doses = list(), inputs = "Cp", known = character(), initial = list(),
    expected = expected(identifiable = c("Kd", "ke0"), non_identifiable = c("Rtot", "ke"), combinations = "Rtot*ke")
  ),
  A10 = list(
    title = "Effect compartment, dynamic receptor binding, linear transduction (known Cp)",
    source = "Janzen et al. 2016, Model 13",
    code = "d/dt(Ce) = ke0*(Cp - Ce)
            d/dt(RC) = kon*(Rtot - RC)*Ce - koff*RC
            effect = ke*RC",
    outputs = "effect", doses = list(), inputs = "Cp", known = character(), initial = list(),
    expected = expected(identifiable = c("ke0", "koff", "kon"), non_identifiable = c("Rtot", "ke"), combinations = "Rtot*ke")
  ),
  A10_Rtot_known = list(
    title = "As A10 with total receptors known (the paper's worked example)",
    source = "Janzen et al. 2016, Model 13, equations 6-8",
    code = "d/dt(Ce) = ke0*(Cp - Ce)
            d/dt(RC) = kon*(Rtot - RC)*Ce - koff*RC
            effect = ke*RC",
    outputs = "effect", doses = list(), inputs = "Cp", known = "Rtot", initial = list(),
    expected = expected(identifiable = c("ke", "ke0", "koff", "kon"))
  )
)
