# Starting points.
#
# Following Algorithm 1 of Sottile and Frumento (2022), a nonparametric estimate
# of F(y | x) is obtained from a piecewise constant hazard model (left
# truncation and right censoring are allowed), and the starting point minimizes
# sum_i w_i (y_i - Q(F(y_i | x_i) | x_i, theta))^2 over 0.1 < F < 0.9.

#' Conditional distribution function from a piecewise-constant hazard model.
#'
#' log h(t | x) = alpha_j + x'beta_j on the interval j (non-proportional hazards,
#' as in package 'pch'), fitted by Poisson regression on split data. If there are
#' too few events, proportional hazards (common beta) are used, and finally a
#' model without covariates.
#' @return F(y_i | x_i).
#' @noRd
pch_cdf <- function(z, y, d, x, w, nbreaks = NULL) {
  n <- length(y)
  x <- as.matrix(x)
  keep <- apply(x, 2, function(v) stats::sd(v) > 0)
  x <- x[, keep, drop = FALSE]
  if (ncol(x)) {
    x <- x[, independent_cols(x), drop = FALSE]
    x <- scale(x)
  }
  p <- ncol(x)
  ev <- y[d == 1]
  if (is.null(nbreaks)) nbreaks <- max(1, min(10, floor(length(ev) / (10 * (p + 1)))))
  br <- unique(stats::quantile(ev, (1:nbreaks) / (nbreaks + 1), names = FALSE))
  br <- c(-Inf, br, Inf)
  k <- length(br) - 1
  zz <- if (all(is.finite(z))) z else rep(-Inf, n)
  # exposure and events by interval
  lo <- pmax(rep(br[-(k + 1)], each = n), zz)
  hi <- pmin(rep(br[-1], each = n), y)
  ref <- min(c(y, zz[is.finite(zz)])) - 1e-8 * diff(range(y))
  lo[!is.finite(lo)] <- ref
  expo <- pmax(hi - lo, 0)
  int <- rep(seq_len(k), each = n)
  ev_j <- as.numeric(d[rep(seq_len(n), k)] == 1 & y[rep(seq_len(n), k)] > rep(br[-(k + 1)], each = n) &
                       y[rep(seq_len(n), k)] <= rep(br[-1], each = n))
  use <- expo > 0
  X <- x[rep(seq_len(n), k), , drop = FALSE]
  I <- outer(int, seq_len(k), "==") + 0
  ww <- w[rep(seq_len(n), k)]
  fitp <- function(M) {
    M <- M[, independent_cols(M[use, , drop = FALSE]), drop = FALSE]
    f <- suppressWarnings(try(stats::glm.fit(M[use, , drop = FALSE], ev_j[use], ww[use],
                                             offset = log(expo[use]), family = stats::poisson(),
                                             control = list(maxit = 50)), silent = TRUE))
    if (inherits(f, "try-error") || !f$converged || any(!is.finite(f$coefficients))) return(NULL)
    b <- f$coefficients; b[is.na(b)] <- 0
    exp(c(M %*% b))
  }
  haz <- NULL
  if (p > 0 && sum(d) >= 10 * k * (p + 1)) {
    M <- cbind(I, do.call(cbind, lapply(seq_len(k), function(j) I[, j] * X)))
    haz <- fitp(M)
  }
  if (is.null(haz) && p > 0) haz <- fitp(cbind(I, X))
  if (is.null(haz)) haz <- fitp(I)
  if (is.null(haz)) {   # Kaplan-Meier (with left truncation)
    km <- km_unif(y, d, w, if (any(is.finite(z))) z else NULL)
    return(pmin(stats::approx(km$u[-1], km$F[-1], xout = y, method = "constant", rule = 2, f = 0)$y, 1 - 1e-3))
  }
  # cumulative hazard from -Inf (or from the first observation) to y
  lo0 <- pmax(rep(br[-(k + 1)], each = n), ref)
  e0 <- pmax(hi - lo0, 0)
  H <- .rowSums(matrix(haz * e0, n, k), n, k)
  1 - exp(-H)
}

# Starting points of Qest for a user-defined Q.
start_generic <- function(z, y, d, x, w, Q, data, start, control) {
  npar <- length(start)
  nostart <- anyNA(start)
  theta <- start; theta[is.na(theta)] <- 0
  if (!nostart && !isTRUE(control$restart)) return(theta)
  Fy <- pch_cdf(z, y, d, x, w)
  use <- which(Fy > 0.1 & Fy < 0.9)
  if (length(use) < npar + 1) use <- seq_along(y)
  du <- data[use, , drop = FALSE]
  obj <- function(th) {
    q <- suppressWarnings(try(callQ(Q, th, Fy[use], du), silent = TRUE))
    if (inherits(q, "try-error") || length(q) != length(use)) return(1e100)
    r <- sum(w[use] * (y[use] - c(q))^2) / length(use) / stats::var(y)
    if (is.finite(r)) r else 1e100
  }
  best <- list(par = theta, value = obj(theta))
  nm <- suppressWarnings(stats::optim(theta, obj, method = "Nelder-Mead",
                                      control = list(maxit = 200 * npar, reltol = 1e-10)))
  if (nm$value < best$value) best <- nm
  bf <- suppressWarnings(try(stats::optim(best$par, obj, method = "BFGS",
                                          control = list(maxit = 200, reltol = 1e-12)), silent = TRUE))
  if (!inherits(bf, "try-error") && bf$value < best$value) best <- bf
  if (best$value >= 1e100) stop("Q(theta, tau, data) could not be evaluated: please provide better starting points",
                                call. = FALSE)
  best$par
}
