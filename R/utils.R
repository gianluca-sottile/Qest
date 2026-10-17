# Small internal utilities.

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x) - log(1 - x)
pmax0 <- function(x) (x + abs(x)) / 2
Ltau <- function(opt, tau) { opt$tau <- tau; opt }

# Cumulative trapezoidal integral of fx on a grid with steps dx (starting from 0).
num.fun <- function(dx, fx) {
  n <- length(dx) + 1
  cumsum(c(0, 0.5 * dx * (fx[1:(n - 1)] + fx[2:n])))
}

# Products X[, j] * X[, h], h >= j (columns ordered as the lower triangle by columns).
tensorX <- function(X) {
  q <- ncol(X)
  out <- NULL
  for (j in seq_len(q)) out <- cbind(out, X[, j] * X[, j:q])
  out
}

# Assemble a symmetric matrix from its lower triangle (by columns).
sym_from_lower <- function(h, npar) {
  H <- matrix(0, npar, npar)
  H[lower.tri(H, diag = TRUE)] <- h
  H <- H + t(H)
  diag(H) <- diag(H) / 2
  H
}

# Linearly independent columns of a model matrix.
independent_cols <- function(x, tol = 1e-6) {
  if (ncol(x) == 0) return(logical(0))
  QR <- qr(x, tol = tol)
  ok <- rep(FALSE, ncol(x))
  ok[QR$pivot[seq_len(QR$rank)]] <- TRUE
  names(ok) <- colnames(x)
  ok
}

#' Frisch-Newton interior point algorithm for linear quantile regression.
#'
#' Pure R implementation of the primal-dual algorithm of Portnoy and Koenker
#' (1997), as in quantreg::rq.fit.fnb. It is only used to compute starting
#' points.
#' @noRd
rq_fnb <- function(x, y, tau = 0.5, beta = 0.99995, eps = 1e-8, maxit = 100) {
  x <- as.matrix(x); n <- nrow(x)
  cc <- -y; b <- (1 - tau) * colSums(x); u <- rep(1, n)
  bound <- function(v, dv) { f <- dv < 0; if (any(f)) min(-v[f] / dv[f]) else 1e20 }
  xx <- rep(1 - tau, n); s <- u - xx
  yy <- qr.coef(qr(x), cc); yy[is.na(yy)] <- 0
  r <- cc - c(x %*% yy)
  z <- pmax(r, 0); w <- z - r
  z <- z + 1e-10 * (z == 0); w <- w + 1e-10 * (w == 0)
  gap <- sum(cc * xx) - sum(yy * b) + sum(w * u)
  it <- 0
  while (gap > eps && it < maxit) {
    it <- it + 1
    q <- 1 / (z / xx + w / s)
    r <- z - w
    R <- chol(crossprod(x * q, x))
    solveR <- function(v) backsolve(R, forwardsolve(t(R), v))
    dy <- solveR(crossprod(x, q * r))
    dx <- q * (c(x %*% dy) - r)
    ds <- -dx
    dz <- -z * (dx / xx + 1)
    dw <- -w * (ds / s + 1)
    fp <- min(beta * min(bound(xx, dx), bound(s, ds)), 1)
    fd <- min(beta * min(bound(w, dw), bound(z, dz)), 1)
    if (min(fp, fd) < 1) {
      mu <- sum(z * xx) + sum(w * s)
      g <- sum((z + fd * dz) * (xx + fp * dx)) + sum((w + fd * dw) * (s + fp * ds))
      mu <- mu * (g / mu)^3 / (2 * n)
      dxdz <- dx * dz; dsdw <- ds * dw
      xinv <- 1 / xx; sinv <- 1 / s
      xi <- mu * (xinv - sinv)
      dy <- solveR(crossprod(x, q * (r + dxdz - dsdw - xi)))
      dx <- q * (c(x %*% dy) + xi - r - dxdz + dsdw)
      ds <- -dx
      dz <- mu * xinv - z - xinv * z * dx - dxdz
      dw <- mu * sinv - w - sinv * w * ds - dsdw
      fp <- min(beta * min(bound(xx, dx), bound(s, ds)), 1)
      fd <- min(beta * min(bound(w, dw), bound(z, dz)), 1)
    }
    xx <- xx + fp * dx; s <- s + fp * ds
    yy <- yy + fd * c(dy); w <- w + fd * dw; z <- z + fd * dz
    gap <- sum(cc * xx) - sum(yy * b) + sum(w * u)
  }
  coef <- -c(yy)
  names(coef) <- colnames(x)
  list(coefficients = coef, residuals = y - c(x %*% coef), n.it = it)
}

# Starting points of linear models (Qnorm, Qunif, Qlm) were computed by the
# Barrodale-Roberts algorithm (Fortran) in Qest <= 1.0.2.
rq.fit.br2 <- function(x, y, tau = 0.5) rq_fnb(x, y, tau)

#' Weights for Truncated Quantile Integration
#'
#' A weighting function \eqn{w(\tau)} that assigns zero weight to the quantiles
#' below `delta.left` and above `1 - delta.right`, to be used as the argument
#' `wtau` of [Qest()], [Qlm()] and [Qcoxph()]. Truncated weights make the
#' estimators more robust to outliers (Sottile and Frumento, 2022, Section 3).
#'
#' @param tau a vector of quantile orders in \eqn{[0, 1]}.
#' @param delta.left,delta.right the proportions of quantiles with zero weight in
#'   the left and right tail.
#' @param smooth logical: if `TRUE`, a smooth version based on the normal
#'   distribution function is used.
#' @param sigma the scale of the smooth version (only used if `smooth = TRUE`).
#'
#' @details
#' With `smooth = FALSE`, \eqn{w(\tau) = I(\delta_L \le \tau \le 1 - \delta_R)}.
#' With `smooth = TRUE`, the indicators are replaced by
#' \eqn{\Phi((\tau - \delta_L)/\sigma)} and \eqn{1 - \Phi((\tau - 1 + \delta_R)/\sigma)}.
#' A `delta` equal to zero means no truncation on that side. The arguments
#' `delta.left`, `delta.right`, `smooth` and `sigma` are passed through the
#' `...` of the fitting functions.
#'
#' @return A vector of weights of the same length as `tau`.
#' @references Sottile, G. and Frumento, P. (2022). Robust estimation and
#'   regression with parametric quantile functions. *Computational Statistics
#'   and Data Analysis*, 171, 107471. \doi{10.1016/j.csda.2022.107471}
#' @examples
#' taus <- seq(0, 1, length.out = 1000)
#' plot(taus, wtrunc(taus, delta.left = 0.05, delta.right = 0.1), type = "l")
#' lines(taus, wtrunc(taus, 0.05, 0.1, smooth = TRUE, sigma = 0.01), col = 2)
#'
#' # in a model
#' set.seed(1)
#' y <- c(rnorm(95), rnorm(5, 10))   # 5% of outliers
#' fit <- Qlm(y ~ 1, wtau = wtrunc, delta.left = 0.05, delta.right = 0.05)
#' coef(fit)
#' @export
wtrunc <- function(tau, delta.left = 0.05, delta.right = 0.05, smooth = FALSE, sigma = 0.01) {
  w1 <- if (delta.right != 0) {
    if (smooth) 1 - stats::pnorm((tau - (1 - delta.right)) / sigma) else as.numeric(tau <= 1 - delta.right)
  } else 1
  w2 <- if (delta.left != 0) {
    if (smooth) stats::pnorm((tau - delta.left) / sigma) else as.numeric(tau >= delta.left)
  } else 1
  w1 * w2
}

callwtau <- function(wtau, tau, opt) {
  opt$tau <- tau
  do.call(wtau, opt)
}

# Split the arguments in ... between wtau and Q.
split_dots <- function(dots, wtau, Q = NULL) {
  argw <- if (!is.null(wtau)) setdiff(names(formals(wtau)), c("tau", "...")) else character(0)
  argQ <- if (is.function(Q)) setdiff(names(formals(Q)), c("theta", "tau", "data", "...")) else character(0)
  nm <- names(dots)
  if (length(dots) && (is.null(nm) || any(nm == "")))
    stop("all arguments in '...' must be named", call. = FALSE)
  unknown <- setdiff(nm, c(argw, argQ))
  if (length(unknown) && !is.null(wtau) && "..." %in% names(formals(wtau))) {
    argw <- c(argw, unknown); unknown <- character(0)
  }
  if (length(unknown))
    stop("unused argument(s) in '...': ", paste(unknown, collapse = ", "),
         " (they are neither arguments of 'wtau' nor of 'Q')", call. = FALSE)
  list(wtau = dots[nm %in% argw], Q = dots[nm %in% argQ])
}
