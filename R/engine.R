# Estimating equations of Qest for a user-defined quantile function.
#
# All integrals in tau are computed on the grid tau$tau (with rectangles on
# [0, tau_1] and [tau_K, 1]). F(y | x) is obtained by inverting the piecewise
# linear interpolant of Q(tau | x) on the grid. The first and second derivatives
# of Q with respect to theta are computed by central differences, and the
# Jacobian is the exact derivative of the discretized estimating equations.
#
# Notation (Sottile and Frumento, 2022, Appendix A):
#   a(tau) = w(tau) Q_theta(tau), A(tau) = int_0^tau a, Ab(tau) = int_0^tau A,
#   U(F) = (Ab(F) - Ab(1)) / (1 - F),
#   g_i = -d_i A(F(y_i)) + (1 - d_i) U(F(y_i)) - U(F(z_i)),
# where F(z_i) = 0, i.e. -U(F(z_i)) = Ab(1), without truncation. With complete
# data, g_i = Ab(1) - A(F(y_i)), the gradient of the loss.

callQ <- function(Q, theta, tau, data) {
  Qoptions <- attr(Q, "Qoptions")
  Qoptions$theta <- theta
  Qoptions$tau <- tau
  if (is.null(attr(Q, "X"))) {
    Qoptions$data <- data
    do.call(Q, Qoptions)
  } else {
    Qoptions$X <- attr(Q, "X")
    do.call(Q, Qoptions)$Q
  }
}

# Q on the grid. Returns NULL if the result is not valid (non finite values,
# decreasing quantiles), which is used to reject a step.
qgrid <- function(Q, theta, TAU, data) {
  q <- suppressWarnings(try(callQ(Q, theta, TAU, data), silent = TRUE))
  if (inherits(q, "try-error")) return(NULL)
  q <- as.matrix(q)
  if (!identical(dim(q), dim(TAU)) || any(!is.finite(q))) return(NULL)
  K <- ncol(q)
  dq <- q[, -1, drop = FALSE] - q[, -K, drop = FALSE]
  if (any(dq < -1e-8 * pmax(1, abs(q[, -1, drop = FALSE])))) return(NULL)
  q
}

# Locate y on the grid: Q[i, k] <= y[i] < Q[i, k + 1].
locate <- function(y, Q0, tg) {
  n <- nrow(Q0); K <- ncol(Q0)
  v <- .rowSums(Q0 <= y, n, K)
  outL <- v == 0; outR <- v == K
  k <- pmin(pmax(v, 1), K - 1)
  i <- seq_len(n)
  QL <- Q0[cbind(i, k)]; QR <- Q0[cbind(i, k + 1)]
  dt <- tg[k + 1] - tg[k]
  dQ <- QR - QL
  s <- ifelse(dQ > 0, (y - QL) / dQ, 0.5)
  s[outL] <- 0; s[outR] <- 1
  s <- pmin(pmax(s, 0), 1)
  slope <- dQ / dt
  out <- outL | outR
  list(k = k, s = s, F = tg[k] + s * dt, dt = dt, slope = slope, out = out, n = n)
}

interp <- function(M, loc) {
  i <- seq_len(loc$n)
  M[cbind(i, loc$k)] * (1 - loc$s) + M[cbind(i, loc$k + 1)] * loc$s
}
islope <- function(M, loc) {
  i <- seq_len(loc$n)
  (M[cbind(i, loc$k + 1)] - M[cbind(i, loc$k)]) / loc$dt
}

# Cumulative integral by row (rectangle on [0, tau_1], trapezoids afterwards).
cumint <- function(a, tg, half_first = FALSE) {
  n <- nrow(a); K <- ncol(a)
  dt <- diff(tg)
  inc <- cbind(tg[1] * a[, 1] * if (half_first) 0.5 else 1,
               (a[, -1, drop = FALSE] + a[, -K, drop = FALSE]) * rep(dt / 2, each = n))
  matrixStats::rowCumsums(inc)
}

# Integrals A, Ab, Ab(1) of a = w * Qtheta.
intblocks <- function(a, tg) {
  A <- cumint(a, tg)
  Ab <- cumint(A, tg, half_first = TRUE)
  K <- length(tg)
  list(A = A, Ab = Ab, Ab1 = Ab[, K] + (1 - tg[K]) * A[, K])
}

# Main function. ctx contains Q, data, y, z, d, w, type, tau, h.
generic_ee <- function(theta, ctx, J = FALSE) {
  tau <- ctx$tau; tg <- tau$tau; TAU <- tau$TAU
  n <- length(ctx$y); npar <- length(theta)
  Q0 <- qgrid(ctx$Q, theta, TAU, ctx$data)
  if (is.null(Q0)) return(NULL)
  h <- ctx$h * pmax(abs(theta), 1)
  Qp <- Qm <- vector("list", npar)
  for (j in seq_len(npar)) {
    e <- replace(numeric(npar), j, h[j])
    Qp[[j]] <- qgrid(ctx$Q, theta + e, TAU, ctx$data)
    Qm[[j]] <- qgrid(ctx$Q, theta - e, TAU, ctx$data)
    if (is.null(Qp[[j]]) || is.null(Qm[[j]])) return(NULL)
  }
  wt <- rep(tau$wtau, each = n)
  ly <- locate(ctx$y, Q0, tg)
  ct <- ctx$type == "ct"
  lz <- if (ct) locate(ctx$z, Q0, tg) else NULL
  d <- ctx$d
  cens <- ctx$type != "u"

  Qth <- B <- vector("list", npar)
  g.i <- matrix(0, n, npar)
  for (j in seq_len(npar)) {
    Qth[[j]] <- (Qp[[j]] - Qm[[j]]) / (2 * h[j])
    B[[j]] <- intblocks(Qth[[j]] * wt, tg)
    Ay <- interp(B[[j]]$A, ly)
    if (!cens) g.i[, j] <- B[[j]]$Ab1 - Ay
    else {
      Uy <- (interp(B[[j]]$Ab, ly) - B[[j]]$Ab1) / (1 - ly$F)
      g.i[, j] <- -d * Ay + (1 - d) * Uy +
        (if (ct) -(interp(B[[j]]$Ab, lz) - B[[j]]$Ab1) / (1 - lz$F) else B[[j]]$Ab1)
    }
  }
  g <- c(ctx$w %*% g.i)
  out <- list(g = g, g.i = g.i, Q0 = Q0, Fy = ly, Fz = lz,
              fy = ifelse(ly$out, 0, 1 / ly$slope), fz = if (ct) ifelse(lz$out, 0, 1 / lz$slope))
  if (!J) return(out)

  # derivatives of F(y) and F(z) w.r.t. theta
  dF <- function(loc) {
    if (is.null(loc)) return(NULL)
    vapply(seq_len(npar), function(k) ifelse(loc$out, 0, -interp(Qth[[k]], loc) / loc$slope), numeric(n))
  }
  Fy_th <- matrix(dF(ly), n); Fz_th <- if (ct) matrix(dF(lz), n)

  # second derivatives of Q
  Q2 <- function(j, k) {
    if (j == k) return((Qp[[j]] - 2 * Q0 + Qm[[j]]) / h[j]^2)
    ej <- replace(numeric(npar), j, h[j]); ek <- replace(numeric(npar), k, h[k])
    qpp <- qgrid(ctx$Q, theta + ej + ek, TAU, ctx$data); qpm <- qgrid(ctx$Q, theta + ej - ek, TAU, ctx$data)
    qmp <- qgrid(ctx$Q, theta - ej + ek, TAU, ctx$data); qmm <- qgrid(ctx$Q, theta - ej - ek, TAU, ctx$data)
    if (is.null(qpp) || is.null(qpm) || is.null(qmp) || is.null(qmm)) return(NULL)
    (qpp - qpm - qmp + qmm) / (4 * h[j] * h[k])
  }
  # derivative of U_j(F) w.r.t. theta_k
  dU <- function(Bj, Bjk, loc, Fk) {
    num <- interp(Bjk$Ab, loc) + islope(Bj$Ab, loc) * Fk - Bjk$Ab1
    num / (1 - loc$F) + (interp(Bj$Ab, loc) - Bj$Ab1) * Fk / (1 - loc$F)^2
  }
  H <- matrix(0, npar, npar)
  for (j in seq_len(npar)) for (k in j:npar) {
    q2 <- Q2(j, k)
    if (is.null(q2)) return(NULL)
    Bjk <- intblocks(q2 * wt, tg)
    pairs <- if (j == k) list(c(j, k)) else list(c(j, k), c(k, j))
    for (jk in pairs) {
      a <- jk[1]; b <- jk[2]   # derivative of g_a w.r.t. theta_b
      dAy <- interp(Bjk$A, ly) + islope(B[[a]]$A, ly) * Fy_th[, b]
      if (!cens) dg <- Bjk$Ab1 - dAy
      else {
        dg <- -d * dAy + (1 - d) * dU(B[[a]], Bjk, ly, Fy_th[, b]) +
          (if (ct) -dU(B[[a]], Bjk, lz, Fz_th[, b]) else Bjk$Ab1)
      }
      H[a, b] <- sum(ctx$w * dg)
    }
  }
  out$J <- H
  out
}

# Loss (the integrated check function), computed on the grid. For censored and
# truncated data, the check function is replaced by its conditional expectation
# (Appendix A), and this is not the function whose gradient is the estimating
# equation.
generic_loss <- function(Q0, y, z, d, w, type, tau, Fy, Fz) {
  tg <- tau$tau; K <- length(tg); dt <- diff(tg)
  tw <- c(tg[1] + dt[1] / 2, (dt[-1] + dt[-(K - 1)]) / 2, dt[K - 1] / 2 + 1 - tg[K]) * tau$wtau
  TAU <- tau$TAU
  dy <- y - Q0
  om <- (dy <= 0) + 0
  if (type != "u") om <- om + (1 - d) * om * (TAU - 1) / (1 - pmin(Fy, 1 - 1e-10))
  if (type == "ct") {
    oz <- (z - Q0 <= 0) + 0
    om <- om - oz - oz * (TAU - 1) / (1 - pmin(Fz, 1 - 1e-10)) + TAU
  }
  sum(w * c((dy * (TAU - om)) %*% tw))
}
