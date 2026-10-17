# Solver of the estimating equations g(theta) = 0, used by all models.
#
# Newton-Raphson with backtracking on ||g||; when the Newton step fails (singular
# Jacobian, or no decrease of ||g||), Levenberg-Marquardt steps
# (J'J + lambda * diag(J'J)) delta = J'g with adaptive lambda. Since
# d/dalpha ||g(theta - alpha J^{-1} g)||^2 = -2 ||g||^2 at alpha = 0, the Newton
# direction is a descent direction for ||g|| whenever J is the Jacobian of g.
#
# eefun(theta, J) must return NULL (invalid theta) or a list with elements g
# (sum over observations), g.i and, if J = TRUE, J (sum over observations).
qest_solve <- function(theta, eefun, n, control, what = "", max.step = 1) {
  tol <- control$tol
  maxit <- control$maxit; if (is.null(maxit)) maxit <- 50 + 25 * length(theta)
  display <- isTRUE(control$display)
  nrm <- function(G) if (is.null(G) || any(!is.finite(G$g))) Inf else sqrt(sum((G$g / n)^2))

  G <- eefun(theta, TRUE)
  if (is.null(G) || any(!is.finite(G$g)) || any(!is.finite(G$J)))
    stop("the estimating equations cannot be evaluated at the starting points", call. = FALSE)
  nm <- nrm(G); nm0 <- nm
  lambda <- 0
  conv <- FALSE
  it <- 0
  if (display) cat(sprintf("%sIter %4d  ||g|| = %.10g\n", what, 0L, nm))
  while (it < maxit) {
    if (nm < tol) { conv <- TRUE; break }
    it <- it + 1
    g <- G$g / n; Jm <- G$J / n
    accepted <- FALSE
    # 1. Newton with backtracking
    # steps are shortened so that |delta_j| <= max.step * max(1, |theta_j|)
    cap <- function(s) { r <- max(abs(s) / (max.step * pmax(1, abs(theta)))); if (r > 1) s / r else s }
    step <- tryCatch(qr.solve(Jm, g, tol = 1e-10), error = function(e) NULL)
    if (!is.null(step) && all(is.finite(step))) {
      step <- cap(step)
      alpha <- 1
      for (b in 1:6) {
        newt <- theta - alpha * step
        G1 <- eefun(newt, FALSE)
        nm1 <- nrm(G1)
        if (nm1 < nm) { accepted <- TRUE; delta <- alpha * step; break }
        alpha <- alpha / 2
      }
    }
    # 2. Levenberg-Marquardt
    if (!accepted) {
      A <- crossprod(Jm); D <- pmax(diag(A), 1e-12 * max(diag(A), 1e-300))
      lambda <- max(lambda, 1e-3)
      while (lambda < 1e12) {
        step <- tryCatch(solve(A + lambda * diag(D, length(D)), crossprod(Jm, g)), error = function(e) NULL)
        if (!is.null(step) && all(is.finite(step))) {
          step <- cap(c(step))
          newt <- theta - c(step)
          G1 <- eefun(newt, FALSE)
          nm1 <- nrm(G1)
          if (nm1 < nm) { accepted <- TRUE; delta <- c(step); lambda <- lambda / 10; break }
        }
        lambda <- lambda * 10
      }
    } else lambda <- 0
    if (!accepted) break
    small <- max(abs(delta) / (abs(theta) + 1)) < 1e-7
    theta <- newt
    G <- eefun(theta, TRUE)
    if (is.null(G) || any(!is.finite(G$J))) stop("the Jacobian cannot be evaluated", call. = FALSE)
    nm <- nrm(G)
    if (display) cat(sprintf("%sIter %4d  ||g|| = %.10g\n", what, it, nm))
    if (small) { conv <- TRUE; break }
  }
  if (nm < tol) conv <- TRUE
  if (!conv && nm < sqrt(tol) * max(1, nm0)) conv <- TRUE   # no further progress at numerical precision
  list(theta = theta, G = G, converged = conv, n.it = it, norm = nm)
}

# Sandwich covariance matrix J^{-1} Omega J^{-T}, Omega = sum_i w_i g_i g_i'.
sandwich <- function(J, g.i, w) {
  Ji <- solve(J)
  Ji %*% crossprod(g.i * sqrt(w)) %*% t(Ji)
}
