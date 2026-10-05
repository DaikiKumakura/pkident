#' Structural identifiability
#'
#' Decides local structural identifiability of the unknown parameters of a
#' [pkpd_model()]. The outputs are differentiated along the model (Lie
#' derivatives), with known input signals contributing their own time
#' derivatives, and evaluated at the known initial state at time zero (the
#' Taylor-series approach). The parameters are locally identifiable if the
#' Jacobian of these output derivatives with respect to the parameters has full
#' rank. The rank is evaluated at several random parameter points on a log
#' scale; a decision is made only if all points agree.
#'
#' Using the known initial state matters for pharmacokinetic models: a bolus
#' dose fixes the initial amount, which is often what makes a volume
#' identifiable. Methods that evaluate the rank at generic initial states ignore
#' this information.
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
#' @param points Number of random parameter points.
#' @param seed Random seed for the points (the global random state is restored).
#' @param max_order Highest derivative order. By default the number of
#'   parameters plus the number of states plus two, capped at 14.
#' @param plateau Number of consecutive orders without a rank increase needed
#'   to conclude rank deficiency.
#' @param tol Relative singular value below which a direction counts as null.
#' @param max_seconds Time limit; beyond it the result is `unresolved`.
#' @param max_exponent,max_support Search limits for monomial combinations.
#' @return An object of class `pkident_structural`. Use [as_table()] to extract
#'   the `"parameters"`, `"combinations"` and `"rank"` tables.
#' @examples
#' r <- structural_identifiability(reference_model("A2"))
#' r
#' as_table(r, "combinations")
#' @export
structural_identifiability <- function(m, points = 5L, seed = 1L, max_order = NULL, plateau = 3L,
                                       tol = 1e-7, max_seconds = 300, max_exponent = 2L, max_support = 4L) {
  if (!inherits(m, "pkident_model")) pki_abort("PKI009", "`m` must be a pkident_model from pkpd_model().")
  started <- Sys.time()
  theta <- m$parameters
  n <- length(theta)
  if (is.null(max_order)) max_order <- min(n + length(m$states) + 2L, 14L)

  # Known input signals: u becomes u_d0; its derivatives are u_d1, u_d2, ...
  in_sym <- function(u, k) paste0("pkiu_", u, "_d", k)
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

  # Random points (log scale for parameters and known constants)
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  pts <- lapply(seq_len(points), function(i) {
    v <- c(stats::setNames(exp(stats::runif(n, log(0.2), log(5))), theta),
           stats::setNames(exp(stats::runif(length(m$known), log(0.5), log(5))), m$known))
    for (u in m$inputs) {
      v[in_sym(u, 0)] <- stats::runif(1, 0.5, 2)
      for (k in seq_len(max_order + 1L)) v[in_sym(u, k)] <- stats::rnorm(1)
    }
    v
  })
  arg_names <- names(pts[[1]])

  h <- outs
  rows <- vector("list", points)
  history <- list()
  status_reason <- NA_character_
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
    ranks <- vapply(seq_len(points), function(i) rank_info(rows[[i]], pts[[i]][theta], tol)$rank, integer(1))
    ambiguous <- vapply(seq_len(points), function(i) rank_info(rows[[i]], pts[[i]][theta], tol)$ambiguous, logical(1))
    history[[length(history) + 1L]] <- data.frame(order = k, point = seq_len(points), rank = ranks, ambiguous = ambiguous)
    final_rank <- ranks
    order_used <- k
    if (all(ranks == n)) break
    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > max_seconds) {
      status_reason <- sprintf("time limit of %s s reached at order %d", max_seconds, k)
      break
    }
  }
  hist <- do.call(rbind, history)

  # Decide
  decision <- "decided"
  if (!is.na(status_reason)) decision <- "unresolved"
  if (any(hist$ambiguous[hist$order == order_used])) {
    decision <- "unresolved"; status_reason <- "singular values fall in the ambiguous band"
  }
  if (length(unique(final_rank)) > 1L) {
    decision <- "unresolved"; status_reason <- "rank differs between random points"
  }
  full <- all(final_rank == n)
  if (decision == "decided" && !full) {
    per_order <- tapply(hist$rank, hist$order, max)
    last <- utils::tail(per_order, plateau + 1L)
    if (length(last) <= plateau || any(diff(last) != 0L)) {
      decision <- "unresolved"
      status_reason <- sprintf("rank still increasing or too few orders (rank %d of %d at order %d)", final_rank[1], n, order_used)
    }
  }

  # Null space and classification
  null_comp <- matrix(NA_real_, points, n, dimnames = list(NULL, theta))
  null_bases <- vector("list", points)
  for (i in seq_len(points)) {
    ri <- rank_info(rows[[i]], pts[[i]][theta], tol)
    null_bases[[i]] <- ri$null
    null_comp[i, ] <- if (ncol(ri$null)) sqrt(rowSums(ri$null^2)) else 0
  }
  comp <- apply(null_comp, 2, max)
  comp_min <- apply(null_comp, 2, min)
  # Adding derivative orders can only shrink the null space, so a parameter
  # outside it at every point is identifiable even if the search stopped early.
  status <- ifelse(comp < 1e-6, "identifiable",
                   ifelse(decision == "decided" & comp_min > 1e-3, "non_identifiable", "unresolved"))

  combos <- data.frame(combination = character(), parameters = character(), stringsAsFactors = FALSE)
  ni <- theta[status == "non_identifiable"]
  if (decision == "decided" && length(ni) >= 2L) {
    combos <- find_combinations(ni, theta, null_bases, max_exponent, min(max_support, length(ni)),
                                n_needed = length(ni) - ncol(null_bases[[1]]))
    in_combo <- unique(unlist(strsplit(combos$parameters, ",", fixed = TRUE)))
    status[theta %in% in_combo & status == "non_identifiable"] <- "combination_only"
  }

  params <- data.frame(parameter = theta, status = unname(status),
                       null_component = unname(comp), stringsAsFactors = FALSE)
  structure(
    list(
      parameters = params,
      combinations = combos,
      rank = hist,
      decision = decision,
      reason = status_reason,
      full_rank = full,
      n_parameters = n,
      order_used = order_used,
      settings = list(method = "lie_taylor", points = points, seed = seed, max_order = max_order,
                      plateau = plateau, tol = tol, max_seconds = max_seconds),
      local = TRUE,
      null_bases = null_bases,
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
  cat("Decision is LOCAL (Lie derivatives at the known initial state; ", x$settings$points,
      " random points, derivative order up to ", x$order_used, ").\n", sep = "")
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
as_table.pkident_structural <- function(x, what = c("parameters", "combinations", "rank"), ...) {
  what <- match.arg(what)
  tibble::as_tibble(x[[what]])
}

# Internals -----------------------------------------------------------------

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

find_combinations <- function(ni, theta, null_bases, max_exponent, max_support, n_needed) {
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
          max(abs(crossprod(N, a))) / sqrt(sum(a^2)) < 1e-6
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
  ok <- TRUE
  for (N in result$null_bases) {
    if (ncol(N) && max(abs(crossprod(N, a))) / sqrt(sum(a^2)) > 1e-6) ok <- FALSE
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
