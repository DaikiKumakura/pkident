# Repeated estimation on simulated data, used to check the expected relative
# standard errors from practical_identifiability(). Also sourced by
# validation/fim-simulation.R.
#
# Data are simulated from the model at the true values with the residual
# error model. Every data set is fitted by maximum likelihood (residual error
# known) on log parameters with Fisher scoring and Levenberg-Marquardt
# damping. All data sets are solved together in one rxode2 call per
# iteration (one subject per data set).

solve_batch <- function(m, model, osf, P, times, input_values) {
  theta <- m$parameters
  n <- length(theta)
  R <- nrow(P)
  vals <- list()
  if (!is.null(model)) {
    sol <- suppressMessages(suppressWarnings(rxode2::rxSolve(
      model, params = as.data.frame(P), events = rxode2::et(times), atol = 1e-10, rtol = 1e-8,
      maxsteps = 1e6, returnType = "data.frame"
    )))
    idcol <- intersect(c("sim.id", "id"), names(sol))
    sid <- if (length(idcol)) sol[[idcol[1]]] else rep(1L, nrow(sol))
    sol <- sol[order(sid, sol$time), ]
    if (nrow(sol) != R * length(times)) stop("solver returned an incomplete grid")
    for (nm in setdiff(osf$args, c(theta, m$known, m$inputs))) vals[[nm]] <- sol[[nm]]
  }
  for (p in c(theta, m$known)) vals[[p]] <- rep(P[, p], each = length(times))
  for (u in m$inputs) vals[[u]] <- rep(input_values[[u]], R)
  res <- matrix(do.call(osf$f, vals[osf$args]), nrow = R * length(times))
  lapply(stats::setNames(seq_along(m$outputs), names(m$outputs)), function(k) {
    first <- (k - 1L) * (n + 1L) + 1L
    list(y = matrix(res[, first], R, byrow = TRUE),                       # R x times
         S = array(res[, first + seq_len(n)], c(length(times), R, n)))   # times x R x n
  })
}

simulate_and_fit <- function(m, design, values, error, inputs = list(), reps = 200, seed = 1, iterations = 60) {
  theta <- m$parameters
  n <- length(theta)
  outs <- names(m$outputs)
  if (is.numeric(design)) design <- stats::setNames(rep(list(design), length(outs)), outs)
  if (inherits(error, "pkident_residual")) error <- stats::setNames(rep(list(error), length(design)), names(design))
  times <- sort(unique(unlist(design)))
  input_lines <- vapply(m$inputs, function(u) paste0(u, " = ", inputs[[u]]), character(1))
  input_values <- stats::setNames(lapply(m$inputs, function(u) {
    rep_len(as.numeric(eval(parse(text = inputs[[u]]), list(t = times), baseenv())), length(times))
  }), m$inputs)
  model <- compile_sensitivity_model(m, input_lines)
  osf <- output_sensitivity_function(m)
  known <- values[m$known]
  params <- function(LT) {
    P <- cbind(exp(LT), matrix(known, nrow(LT), length(known), byrow = TRUE))
    colnames(P) <- c(theta, m$known)
    P
  }

  set.seed(seed)
  truth <- solve_batch(m, model, osf, params(matrix(log(values[theta]), 1)), times, input_values)
  obs <- lapply(names(design), function(o) {
    idx <- match(design[[o]], times)
    f <- truth[[o]]$y[1, idx]
    sdv <- sqrt(error[[o]]$add^2 + (error[[o]]$prop * f)^2)
    matrix(f, reps, length(idx), byrow = TRUE) + matrix(stats::rnorm(reps * length(idx)), reps) *
      matrix(sdv, reps, length(idx), byrow = TRUE)
  })
  names(obs) <- names(design)

  # -2 log-likelihood, score and expected information for every data set
  evaluate <- function(LT) {
    sol <- solve_batch(m, model, osf, params(LT), times, input_values)
    obj <- numeric(reps)
    score <- matrix(0, reps, n)
    info <- array(0, c(n, n, reps))
    for (o in names(design)) {
      idx <- match(design[[o]], times)
      e <- error[[o]]
      f <- sol[[o]]$y[, idx, drop = FALSE]
      v <- e$add^2 + (e$prop * f)^2
      r <- obs[[o]] - f
      obj <- obj + rowSums(log(v) + r^2 / v)
      for (i in seq_len(reps)) {
        S <- sol[[o]]$S[idx, i, , drop = TRUE]
        S <- matrix(S, length(idx), n) * matrix(exp(LT[i, ]), length(idx), n, byrow = TRUE)   # log scale
        dv <- 2 * e$prop^2 * f[i, ]
        score[i, ] <- score[i, ] + colSums(S * ((1 / v[i, ] - r[i, ]^2 / v[i, ]^2) * dv - 2 * r[i, ] / v[i, ]))
        w <- 1 / v[i, ] + 2 * e$prop^4 * f[i, ]^2 / v[i, ]^2
        info[, , i] <- info[, , i] + 2 * crossprod(S * sqrt(w))   # expected Hessian of -2LL
      }
    }
    obj[!is.finite(obj)] <- Inf
    list(obj = obj, score = score, info = info)
  }

  LT <- matrix(log(values[theta]), reps, n, byrow = TRUE) + matrix(stats::rnorm(reps * n, 0, 0.1), reps)
  lambda <- rep(1e-3, reps)
  cur <- evaluate(LT)
  done <- rep(FALSE, reps)
  for (it in seq_len(iterations)) {
    step <- matrix(0, reps, n)
    for (i in which(!done)) {
      H <- cur$info[, , i]
      step[i, ] <- tryCatch(-solve(H + lambda[i] * diag(diag(H), n), cur$score[i, ]), error = function(e) rep(0, n))
    }
    step <- pmax(pmin(step, 1), -1)
    trial <- evaluate(LT + step)
    better <- trial$obj < cur$obj & !done
    small <- abs(cur$obj - trial$obj) < 1e-10 | rowSums(abs(step)) < 1e-9
    LT[better, ] <- LT[better, ] + step[better, ]
    cur$obj[better] <- trial$obj[better]
    cur$score[better, ] <- trial$score[better, ]
    cur$info[, , better] <- trial$info[, , better]
    lambda[better] <- lambda[better] / 10
    lambda[!better] <- lambda[!better] * 10
    done <- done | (small & !better) | (better & small)
    if (all(done)) break
  }
  # converged: gradient small relative to the information
  converged <- vapply(seq_len(reps), function(i) {
    g <- cur$score[i, ]
    H <- cur$info[, , i]
    q <- tryCatch(sum(g * solve(H, g)), error = function(e) Inf)
    is.finite(cur$obj[i]) && q < 1e-6
  }, logical(1))
  est <- exp(LT)
  colnames(est) <- theta
  list(estimates = est, converged = converged, iterations = it)
}

# Empirical SD of log estimates (the relative SE) next to the expected RSE.
compare_rse <- function(p, sim) {
  le <- log(sim$estimates[sim$converged, , drop = FALSE])
  emp <- 100 * apply(le, 2, stats::sd)
  data.frame(parameter = p$parameters$parameter, expected_rse = p$parameters$rse_percent,
             empirical_rse = unname(emp[p$parameters$parameter]),
             ratio = unname(emp[p$parameters$parameter]) / p$parameters$rse_percent,
             converged = sum(sim$converged), reps = length(sim$converged))
}
