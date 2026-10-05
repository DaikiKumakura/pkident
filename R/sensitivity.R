# Numerical check from output sensitivities (milestone S4).
#
# The forward sensitivity equations dS/dt = (df/dx) S + df/dtheta, with
# S(0) = dx(0)/dtheta, are generated symbolically and solved with rxode2 at the
# same random parameter points as the Lie method. The output sensitivities on a
# dense time grid form a matrix whose numerical rank is compared with the rank
# of the Lie-derivative Jacobian. Known input signals are replaced by a smooth
# random signal (a constant plus three sinusoids).

sens_state <- function(s, p) paste0("pkis_", s, "__", p)
input_coef <- function(u, what, j = "") paste0("pkif_", u, "_", what, j)

# A random smooth positive input signal for each known input.
draw_input_signals <- function(inputs, points) {
  lapply(seq_len(points), function(i) {
    v <- numeric()
    for (u in inputs) {
      v[input_coef(u, "b")] <- stats::runif(1, 1, 2)
      for (j in 1:3) {
        v[input_coef(u, "a", j)] <- stats::runif(1, 0.1, 0.3)
        v[input_coef(u, "w", j)] <- exp(stats::runif(1, log(0.2), log(3)))
        v[input_coef(u, "p", j)] <- stats::runif(1, 0, 2 * pi)
      }
    }
    v
  })
}

input_signal_code <- function(u) {
  terms <- vapply(1:3, function(j) {
    sprintf("%s*sin(%s*t + %s)", input_coef(u, "a", j), input_coef(u, "w", j), input_coef(u, "p", j))
  }, character(1))
  paste0(u, " = ", input_coef(u, "b"), " + ", paste(terms, collapse = " + "))
}

input_signal_value <- function(u, v, times) {
  out <- v[[input_coef(u, "b")]]
  for (j in 1:3) {
    out <- out + v[[input_coef(u, "a", j)]] * sin(v[[input_coef(u, "w", j)]] * times + v[[input_coef(u, "p", j)]])
  }
  out
}

# rxode2 code for the states and their parameter sensitivities.
sensitivity_code <- function(m, input_lines = vapply(m$inputs, input_signal_code, character(1))) {
  theta <- m$parameters
  lines <- unname(input_lines)
  for (s in m$states) {
    lines <- c(lines, sprintf("d/dt(%s) = %s", s, as.character(m$odes[[s]])))
    if (as.character(m$initial[[s]]) != "0") {
      lines <- c(lines, sprintf("%s(0) = %s", s, as.character(m$initial[[s]])))
    }
  }
  for (s in m$states) {
    f <- m$odes[[s]]
    for (p in theta) {
      rhs <- symengine::D(f, p)
      for (k in m$states) rhs <- rhs + symengine::D(f, k) * symengine::S(sens_state(k, p))
      lines <- c(lines, sprintf("d/dt(%s) = %s", sens_state(s, p), as.character(symengine::expand(rhs))))
      s0 <- symengine::D(m$initial[[s]], p)
      if (as.character(s0) != "0") lines <- c(lines, sprintf("%s(0) = %s", sens_state(s, p), as.character(s0)))
    }
  }
  paste(lines, collapse = "\n")
}

# Output sensitivities dy/dtheta as functions of states, sensitivities,
# parameters, known constants and input values (vectorised over time).
output_sensitivity_function <- function(m) {
  theta <- m$parameters
  exprs <- list()
  for (o in names(m$outputs)) {
    h <- m$outputs[[o]]
    exprs[[length(exprs) + 1L]] <- h
    for (p in theta) {
      e <- symengine::D(h, p)
      for (k in m$states) e <- e + symengine::D(h, k) * symengine::S(sens_state(k, p))
      exprs[[length(exprs) + 1L]] <- e
    }
  }
  args <- c(m$states, unlist(lapply(m$states, function(s) sens_state(s, theta))), theta, m$known, m$inputs)
  list(f = as.function(do.call(symengine::Vector, exprs), args = args), args = args)
}

compile_sensitivity_model <- function(m, input_lines = vapply(m$inputs, input_signal_code, character(1))) {
  if (!length(m$states)) return(NULL)
  tryCatch(
    suppressMessages(rxode2::rxode2(sensitivity_code(m, input_lines))),
    error = function(e) pki_abort("PKI010", paste("rxode2 could not compile the sensitivity equations:", conditionMessage(e)))
  )
}

# Outputs and their sensitivities dy/dtheta (natural scale) at the given
# times. `v` holds parameter values, known constants and any input-signal
# coefficients; `input_values` is a named list of input values at `times`.
solve_output_sensitivities <- function(m, model, osf, v, times, input_values, rtol, atol) {
  theta <- m$parameters
  n <- length(theta)
  vals <- list()
  if (!is.null(model)) {
    sol <- tryCatch(
      suppressMessages(suppressWarnings(rxode2::rxSolve(
        model, params = v, events = rxode2::et(times), atol = atol, rtol = rtol,
        maxsteps = 1e6, returnType = "data.frame"
      ))),
      error = function(e) pki_abort("PKI010", paste("rxode2 could not solve the sensitivity equations:", conditionMessage(e)))
    )
    if (nrow(sol) != length(times)) pki_abort("PKI010", "The sensitivity equations could not be solved on the time grid.")
    for (nm in setdiff(osf$args, c(theta, m$known, m$inputs))) vals[[nm]] <- sol[[nm]]
  }
  for (p in c(theta, m$known)) vals[[p]] <- rep(v[[p]], length(times))
  for (u in m$inputs) vals[[u]] <- input_values[[u]]
  res <- matrix(do.call(osf$f, vals[osf$args]), nrow = length(times))
  if (any(!is.finite(res))) pki_abort("PKI010", "Non-finite outputs or output sensitivities.")
  out <- lapply(seq_along(m$outputs), function(k) {
    first <- (k - 1L) * (n + 1L) + 1L
    list(y = res[, first], S = res[, first + seq_len(n), drop = FALSE])
  })
  stats::setNames(out, names(m$outputs))
}

# Rank of the scaled output sensitivity matrix at each point.
sensitivity_ranks <- function(m, pts, signals, times, tol, rtol, atol) {
  theta <- m$parameters
  n <- length(theta)
  osf <- output_sensitivity_function(m)
  model <- compile_sensitivity_model(m)
  lapply(seq_along(pts), function(i) {
    v <- c(pts[[i]][c(theta, m$known)], signals[[i]])
    inp <- stats::setNames(lapply(m$inputs, input_signal_value, v = signals[[i]], times = times), m$inputs)
    sol <- solve_output_sensitivities(m, model, osf, v, times, inp, rtol, atol)
    blocks <- lapply(sol, function(o) {
      Sy <- sweep(o$S, 2, pts[[i]][theta], `*`)   # log-parameter scale
      scale <- max(abs(o$y), abs(Sy))
      if (!is.finite(scale) || scale == 0) scale <- 1
      Sy / scale
    })
    J <- do.call(rbind, blocks)
    sv <- svd(J, nu = 0, nv = n)
    d <- c(sv$d, rep(0, max(0, n - length(sv$d))))
    ratio <- d / max(d[1], .Machine$double.xmin)
    r <- sum(ratio > tol)
    list(rank = as.integer(r), ratio = ratio,
         ambiguous = any(ratio <= tol & ratio > tol * 1e-3),
         null = sv$v[, seq_len(n) > r, drop = FALSE])
  })
}

# Largest principal-angle sine between two subspaces (0 = identical).
subspace_distance <- function(A, B) {
  if (ncol(A) != ncol(B)) return(1)
  if (!ncol(A)) return(0)
  qa <- qr.Q(qr(A)); qb <- qr.Q(qr(B))
  s <- svd(crossprod(qa, qb), nu = 0, nv = 0)$d
  sqrt(max(0, 1 - min(s)^2))
}
