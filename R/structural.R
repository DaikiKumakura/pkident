#' Structural identifiability
#'
#' Decides local structural identifiability of the unknown parameters of a
#' [pkpd_model()] with two independent methods.
#'
#' **Lie method (primary).** The outputs are differentiated along the model
#' (Lie derivatives), with known input signals contributing their own time
#' derivatives, and evaluated at the known initial state at time zero (the
#' Taylor-series approach). The parameters are locally identifiable if the
#' Jacobian of these output derivatives with respect to the parameters has full
#' rank. Using the known initial state matters for pharmacokinetic models: a
#' bolus dose fixes the initial amount, which is often what makes a volume
#' identifiable. Methods that evaluate the rank at generic initial states ignore
#' this information.
#'
#' **Sensitivity method (numerical check).** The forward sensitivity equations
#' are generated symbolically and solved with rxode2 at the same parameter
#' points; known input signals are replaced by a smooth random signal. The rank
#' of the noise-free output sensitivity matrix on a dense time grid is compared
#' with the Lie rank, and the null spaces of the two methods are compared at
#' every point.
#'
#' With `method = "both"` (the default), a decision is made only if the two
#' methods agree at every point; otherwise the result is `unresolved`. Ranks
#' are evaluated at several random parameter points on a log scale, and a
#' decision is made only if all points agree.
#'
#' Rank-deficient parameters are reported with the directions along which they
#' can change together without changing the outputs. Where these directions
#' correspond to products of powers of parameters (for example `V/F`), the
#' identifiable combinations are reported.
#'
#' Decisions are **local**: a parameter reported as identifiable may still have
#' a finite number of alternative values elsewhere in parameter space (for
#' example the flip-flop of absorption and elimination rates).
#'
#' @param m A `pkident_model` from [pkpd_model()].
#' @param method `"both"` (Lie method checked by sensitivities), `"lie"` or
#'   `"sensitivity"`.
#' @param points Number of random parameter points.
#' @param seed Random seed for the points (the global random state is restored).
#' @param max_order Highest derivative order of the Lie method. By default the
#'   number of parameters plus the number of states plus two, capped at 14.
#' @param plateau Number of consecutive orders without a rank increase needed
#'   to conclude rank deficiency (Lie method).
#' @param tol Relative singular value below which a direction counts as null
#'   (Lie method).
#' @param max_seconds Time limit of the Lie method; beyond it the result is
#'   `unresolved`.
#' @param max_exponent,max_support Search limits for monomial combinations.
#' @param times Time grid of the sensitivity method. By default 0 and 150
#'   log-spaced times from 0.001 to 200, which covers the rate constants
#'   implied by the random parameter points (between 0.2 and 5).
#' @param sens_tol Relative singular value below which a direction counts as
#'   null (sensitivity method).
#' @param rtol,atol Solver tolerances of the sensitivity method.
#' @return An object of class `pkident_structural`. Use [as_table()] to extract
#'   the `"parameters"`, `"combinations"`, `"rank"` (Lie method) and
#'   `"agreement"` (comparison of the methods at each point) tables.
#' @examples
#' r <- structural_identifiability(reference_model("A2"), method = "lie")
#' r
#' as_table(r, "combinations")
#' @export
structural_identifiability <- function(m, method = c("both", "lie", "sensitivity"), points = 5L, seed = 1L,
                                       max_order = NULL, plateau = 3L, tol = 1e-7, max_seconds = 300,
                                       max_exponent = 2L, max_support = 4L, times = NULL,
                                       sens_tol = 1e-7, rtol = 1e-10, atol = 1e-12) {
  if (!inherits(m, "pkident_model")) pki_abort("PKI009", "`m` must be a pkident_model from pkpd_model().")
  method <- match.arg(method)
  started <- Sys.time()
  theta <- m$parameters
  n <- length(theta)
  if (is.null(max_order)) max_order <- min(n + length(m$states) + 2L, 14L)
  if (is.null(times)) times <- c(0, exp(seq(log(1e-3), log(200), length.out = 150L)))
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) || any(times < 0) || is.unsorted(times, strictly = TRUE)) {
    pki_abort("PKI009", "`times` must be increasing, finite and non-negative.")
  }
  in_sym <- function(u, k) paste0("pkiu_", u, "_d", k)

  # Random points (log scale for parameters and known constants)
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  pts <- draw_points(m, points, max_order, in_sym)
  signals <- draw_input_signals(m$inputs, points)

  lie <- NULL
  sens <- NULL
  if (method != "sensitivity") {
    lie <- lie_analysis(m, pts, in_sym, max_order, plateau, tol, max_seconds, started)
    lie$status <- classify_parameters(theta, lie$null_bases, lie$decision, id_tol = 1e-6, ni_tol = 1e-3)
  }
  if (method != "lie") {
    sens <- sensitivity_analysis(m, pts, signals, times, sens_tol, rtol, atol)
    sens$status <- classify_parameters(theta, sens$null_bases, sens$decision, id_tol = 1e-4, ni_tol = 1e-2)
  }

  # Reconcile
  agreement <- NULL
  if (method == "both") {
    dist <- vapply(seq_len(points), function(i) subspace_distance(lie$null_bases[[i]], sens$null_bases[[i]]), numeric(1))
    agreement <- data.frame(point = seq_len(points), lie_rank = lie$ranks, sensitivity_rank = sens$ranks,
                            null_space_distance = dist, agree = lie$ranks == sens$ranks & dist < 1e-4)
    primary <- lie
    decision <- if (lie$decision == "decided" && sens$decision == "decided" && all(agreement$agree)) "decided" else "unresolved"
    reasons <- c(lie$reason, sens$reason,
                 if (!all(agreement$agree)) sprintf("Lie and sensitivity methods disagree at point(s) %s",
                                                    paste(which(!agreement$agree), collapse = ", ")))
    reason <- if (decision == "decided") NA_character_ else paste(stats::na.omit(reasons), collapse = "; ")
    status <- ifelse(lie$status == sens$status, lie$status, "unresolved")
    ctol <- 1e-6
  } else {
    primary <- if (method == "lie") lie else sens
    decision <- primary$decision
    reason <- primary$reason
    status <- primary$status
    ctol <- if (method == "lie") 1e-6 else 1e-4
  }

  combos <- data.frame(combination = character(), parameters = character(), stringsAsFactors = FALSE)
  ni <- theta[status == "non_identifiable"]
  if (decision == "decided" && length(ni) >= 2L) {
    combos <- find_combinations(ni, theta, primary$null_bases, max_exponent, min(max_support, length(ni)),
                                n_needed = length(ni) - ncol(primary$null_bases[[1]]), ctol = ctol)
    in_combo <- unique(unlist(strsplit(combos$parameters, ",", fixed = TRUE)))
    status[theta %in% in_combo & status == "non_identifiable"] <- "combination_only"
  }

  params <- data.frame(parameter = theta, status = unname(status),
                       null_component = unname(primary$null_component), stringsAsFactors = FALSE)
  structure(
    list(
      parameters = params,
      combinations = combos,
      rank = if (is.null(lie)) data.frame(order = integer(), point = integer(), rank = integer(), ambiguous = logical()) else lie$history,
      agreement = agreement,
      decision = decision,
      reason = reason,
      full_rank = all(primary$ranks == n),
      n_parameters = n,
      order_used = if (is.null(lie)) NA_integer_ else lie$order_used,
      settings = list(method = method, points = points, seed = seed, max_order = max_order,
                      plateau = plateau, tol = tol, max_seconds = max_seconds, sens_tol = sens_tol,
                      rtol = rtol, atol = atol, n_times = length(times), combination_tol = ctol),
      local = TRUE,
      null_bases = primary$null_bases,
      sensitivity_ratios = if (is.null(sens)) NULL else sens$ratios,
      model = m,
      versions = list(pkident = as.character(utils::packageVersion("pkident")),
                      rxode2 = m$source$rxode2_version, symengine = m$source$symengine_version),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ),
    class = "pkident_structural"
  )
}

#' @export
print.pkident_structural <- function(x, ...) {
  cat("<pkident structural identifiability>\n")
  s <- x$settings
  how <- switch(s$method,
    both = sprintf("Lie derivatives at the known initial state up to order %d, checked by output sensitivities", x$order_used),
    lie = sprintf("Lie derivatives at the known initial state up to order %d", x$order_used),
    sensitivity = sprintf("output sensitivities on %d time points", s$n_times))
  cat("Decision is LOCAL (", how, "; ", s$points, " random points).\n", sep = "")
  if (!is.null(x$agreement)) {
    cat(sprintf("Methods agree at %d of %d points.\n", sum(x$agreement$agree), nrow(x$agreement)))
  }
  if (x$decision == "unresolved") cat("Unresolved: ", x$reason, "\n", sep = "")
  for (i in seq_len(nrow(x$parameters))) {
    cat(sprintf("  %-12s %s\n", x$parameters$parameter[i], x$parameters$status[i]))
  }
  if (nrow(x$combinations)) {
    cat("Identifiable combinations:\n")
    for (cmb in x$combinations$combination) cat("  ", cmb, "\n", sep = "")
  }
  invisible(x)
}

#' Extract result tables
#'
#' @param x A pkident result.
#' @param what Which table.
#' @param ... Unused.
#' @return A tibble.
#' @export
as_table <- function(x, what, ...) UseMethod("as_table")

#' @rdname as_table
#' @export
as_table.pkident_structural <- function(x, what = c("parameters", "combinations", "rank", "agreement"), ...) {
  what <- match.arg(what)
  if (what == "agreement" && is.null(x$agreement)) {
    pki_abort("PKI009", "The agreement table needs method = \"both\".")
  }
  tibble::as_tibble(x[[what]])
}

# Lie method ----------------------------------------------------------------

lie_analysis <- function(m, pts, in_sym, max_order, plateau, tol, max_seconds, started) {
  theta <- m$parameters
  n <- length(theta)
  points <- length(pts)

  # Known input signals: u becomes u_d0; its derivatives are u_d1, u_d2, ...
  to_input <- function(e) {
    for (u in m$inputs) e <- symengine::subs(e, symengine::S(u), symengine::S(in_sym(u, 0)))
    e
  }
  odes <- lapply(m$odes, to_input)
  outs <- lapply(m$outputs, to_input)
  init <- lapply(m$initial, to_input)
  lie <- function(h, k) {
    out <- symengine::S("0")
    for (s in m$states) out <- out + symengine::D(h, s) * odes[[s]]
    for (u in m$inputs) for (j in 0:k) out <- out + symengine::D(h, in_sym(u, j)) * symengine::S(in_sym(u, j + 1))
    symengine::expand(out)
  }
  at_initial <- function(e) {
    for (s in m$states) e <- symengine::subs(e, symengine::S(s), init[[s]])
    e
  }
  arg_names <- names(pts[[1]])

  h <- outs
  rows <- vector("list", points)
  history <- list()
  reason <- NA_character_
  final_rank <- rep(NA_integer_, points)
  order_used <- NA_integer_
  for (k in 0:max_order) {
    if (k > 0) h <- lapply(h, lie, k = k - 1L)
    exprs <- list()
    for (o in names(h)) {
      Hk <- at_initial(h[[o]])
      for (p in theta) exprs[[length(exprs) + 1L]] <- symengine::D(Hk, p)
    }
    f <- as.function(do.call(symengine::Vector, exprs), args = arg_names)
    for (i in seq_len(points)) {
      val <- do.call(f, as.list(pts[[i]][arg_names]))
      block <- matrix(as.numeric(val), nrow = length(h), byrow = TRUE)
      rows[[i]] <- rbind(rows[[i]], block)
    }
    info <- lapply(seq_len(points), function(i) rank_info(rows[[i]], pts[[i]][theta], tol))
    ranks <- vapply(info, `[[`, integer(1), "rank")
    ambiguous <- vapply(info, `[[`, logical(1), "ambiguous")
    history[[length(history) + 1L]] <- data.frame(order = k, point = seq_len(points), rank = ranks, ambiguous = ambiguous)
    final_rank <- ranks
    order_used <- k
    if (all(ranks == n)) break
    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > max_seconds) {
      reason <- sprintf("time limit of %s s reached at order %d", max_seconds, k)
      break
    }
  }
  hist <- do.call(rbind, history)

  decision <- "decided"
  if (!is.na(reason)) decision <- "unresolved"
  if (any(hist$ambiguous[hist$order == order_used])) {
    decision <- "unresolved"; reason <- "singular values fall in the ambiguous band"
  }
  if (length(unique(final_rank)) > 1L) {
    decision <- "unresolved"; reason <- "rank differs between random points"
  }
  if (decision == "decided" && !all(final_rank == n)) {
    per_order <- tapply(hist$rank, hist$order, max)
    last <- utils::tail(per_order, plateau + 1L)
    if (length(last) <= plateau || any(diff(last) != 0L)) {
      decision <- "unresolved"
      reason <- sprintf("rank still increasing or too few orders (rank %d of %d at order %d)", final_rank[1], n, order_used)
    }
  }
  null_bases <- lapply(info, `[[`, "null")
  list(ranks = final_rank, null_bases = null_bases, null_component = null_component(null_bases, theta),
       decision = decision, reason = reason, history = hist, order_used = order_used)
}

# Sensitivity method --------------------------------------------------------

sensitivity_analysis <- function(m, pts, signals, times, sens_tol, rtol, atol) {
  info <- sensitivity_ranks(m, pts, signals, times, sens_tol, rtol, atol)
  ranks <- vapply(info, `[[`, integer(1), "rank")
  decision <- "decided"
  reason <- NA_character_
  if (any(vapply(info, `[[`, logical(1), "ambiguous"))) {
    decision <- "unresolved"; reason <- "sensitivity singular values fall in the ambiguous band"
  }
  if (length(unique(ranks)) > 1L) {
    decision <- "unresolved"; reason <- "sensitivity rank differs between random points"
  }
  null_bases <- lapply(info, `[[`, "null")
  list(ranks = ranks, null_bases = null_bases, null_component = null_component(null_bases, m$parameters),
       decision = decision, reason = reason, ratios = lapply(info, `[[`, "ratio"))
}

# Shared --------------------------------------------------------------------

null_component <- function(null_bases, theta) {
  comp <- vapply(null_bases, function(N) if (ncol(N)) sqrt(rowSums(N^2)) else rep(0, length(theta)),
                 numeric(length(theta)))
  apply(matrix(comp, nrow = length(theta)), 1, max)
}

# A parameter outside the null space at every point is identifiable (adding
# derivative orders or time points can only shrink the null space); one inside
# it at every point is not identifiable alone, provided a decision was made.
classify_parameters <- function(theta, null_bases, decision, id_tol, ni_tol) {
  comp <- vapply(null_bases, function(N) if (ncol(N)) sqrt(rowSums(N^2)) else rep(0, length(theta)),
                 numeric(length(theta)))
  comp <- matrix(comp, nrow = length(theta))
  cmax <- apply(comp, 1, max)
  cmin <- apply(comp, 1, min)
  ifelse(cmax < id_tol, "identifiable",
         ifelse(decision == "decided" & cmin > ni_tol, "non_identifiable", "unresolved"))
}

# Internals -----------------------------------------------------------------

draw_points <- function(m, points, max_order, in_sym) {
  n <- length(m$parameters)
  lapply(seq_len(points), function(i) {
    v <- c(stats::setNames(exp(stats::runif(n, log(0.2), log(5))), m$parameters),
           stats::setNames(exp(stats::runif(length(m$known), log(0.5), log(5))), m$known))
    for (u in m$inputs) {
      v[in_sym(u, 0)] <- stats::runif(1, 0.5, 2)
      for (k in seq_len(max_order + 1L)) v[in_sym(u, k)] <- stats::rnorm(1)
    }
    v
  })
}

rank_info <- function(J, theta_values, tol) {
  Js <- sweep(J, 2, theta_values, `*`)             # log-parameter scale
  nrm <- apply(abs(Js), 1, max)
  Js <- Js[nrm > 0, , drop = FALSE]
  n <- ncol(J)
  if (!nrow(Js)) return(list(rank = 0L, ambiguous = FALSE, null = diag(n), ratio = numeric()))
  Js <- Js / apply(abs(Js), 1, max)
  sv <- svd(Js, nu = 0, nv = n)
  d <- c(sv$d, rep(0, max(0, n - length(sv$d))))
  ratio <- d / max(d[1], .Machine$double.xmin)
  r <- sum(ratio > tol)
  ambiguous <- any(ratio <= tol & ratio > tol * 1e-4)
  list(rank = as.integer(r), ambiguous = ambiguous, null = sv$v[, seq_len(n) > r, drop = FALSE], ratio = ratio)
}

find_combinations <- function(ni, theta, null_bases, max_exponent, max_support, n_needed, ctol = 1e-6) {
  idx <- match(ni, theta)
  cands <- list()
  for (sz in 2:max_support) {
    for (sub in utils::combn(length(ni), sz, simplify = FALSE)) {
      grid <- as.matrix(expand.grid(rep(list(c(-max_exponent:-1, 1:max_exponent)), sz)))
      for (g in seq_len(nrow(grid))) {
        a <- numeric(length(theta)); a[idx[sub]] <- grid[g, ]
        if (a[a != 0][1] < 0) next                      # canonical sign
        ok <- all(vapply(null_bases, function(N) {
          if (!ncol(N)) return(TRUE)
          max(abs(crossprod(N, a))) / sqrt(sum(a^2)) < ctol
        }, logical(1)))
        if (ok) cands[[length(cands) + 1L]] <- a
      }
    }
  }
  chosen <- list()
  if (length(cands)) {
    key <- vapply(cands, function(a) sum(a != 0) * 100 + sum(abs(a)), numeric(1))
    for (a in cands[order(key)]) {
      M <- do.call(rbind, c(chosen, list(a)))
      if (qr(M)$rank > length(chosen)) chosen[[length(chosen) + 1L]] <- a
      if (length(chosen) >= n_needed) break
    }
  }
  data.frame(
    combination = vapply(chosen, monomial_string, character(1), theta = theta),
    parameters = vapply(chosen, function(a) paste(theta[a != 0], collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
}

monomial_string <- function(a, theta) {
  term <- function(p, e) if (e == 1) p else paste0(p, "^", e)
  num <- vapply(which(a > 0), function(j) term(theta[j], a[j]), character(1))
  den <- vapply(which(a < 0), function(j) term(theta[j], -a[j]), character(1))
  top <- if (length(num)) paste(num, collapse = "*") else "1"
  if (!length(den)) return(top)
  paste0(top, "/", if (length(den) > 1) paste0("(", paste(den, collapse = "*"), ")") else den)
}

# Is a monomial combination (for example "V/F") identifiable under a result?
combination_identifiable <- function(result, combination) {
  theta <- result$parameters$parameter
  a <- parse_monomial(combination, theta)
  ctol <- result$settings$combination_tol
  ok <- TRUE
  for (N in result$null_bases) {
    if (ncol(N) && max(abs(crossprod(N, a))) / sqrt(sum(a^2)) > ctol) ok <- FALSE
  }
  ok
}

parse_monomial <- function(s, theta) {
  a <- stats::setNames(numeric(length(theta)), theta)
  parts <- strsplit(gsub("[()]", "", s), "/", fixed = TRUE)[[1]]
  add <- function(txt, sign) {
    for (f in strsplit(txt, "*", fixed = TRUE)[[1]]) {
      f <- trimws(f); if (f == "1" || !nzchar(f)) next
      pe <- strsplit(f, "^", fixed = TRUE)[[1]]
      if (!pe[1] %in% theta) pki_abort("PKI009", sprintf("`%s` is not a parameter.", pe[1]))
      a[pe[1]] <<- a[pe[1]] + sign * (if (length(pe) > 1) as.numeric(pe[2]) else 1)
    }
  }
  add(parts[1], 1)
  if (length(parts) > 1) add(parts[2], -1)
  unname(a)
}
