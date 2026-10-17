# Continuous Poisson family (Ilienko, 2013) for jittered counts: log(lambda) = x'theta.

#### Poisson family

tau.pois <- function(tau){

  tau$A_0.1_1 <- matrix(
    scan(text = "2.001821e+00  3.942654e-02  1.024726e+00  6.644414e-02 -1.044504e+00  1.405770e-02
              -9.507657e-06 -2.932617e-01  3.306340e-01 -1.018398e+00  1.089903e-02  9.311177e-01
              -3.244777e-02  3.795771e-05  9.876184e-01 -6.101660e-01  1.081108e+00  2.056688e-02
              -1.139983e+00  3.128521e-02 -3.465038e-05 -8.387406e-03  4.954290e-02  2.501427e-02
              -5.755154e-03  1.665659e-03  1.426357e-03 -2.048700e-06 -4.646329e-02  4.284785e-02
              -6.034803e-02 -2.021200e-03  6.742223e-02 -1.645966e-03  1.762123e-06 -1.216432e-03
               1.659728e-03 -7.779998e-04 -1.613607e-04  1.399646e-03 -2.580111e-06 -9.005382e-09
              -3.341377e-07  4.212620e-07 -2.238117e-07 -4.083766e-08  3.785965e-07 -1.428820e-09
              -1.379904e-12 -4.824025e-04  4.797075e-04 -5.027079e-04 -3.179307e-05  6.249717e-04
              -1.193474e-05  1.166221e-08 -1.184018e-07  1.187020e-07 -1.178758e-07 -8.238543e-09
               1.497869e-07 -2.701888e-09  2.571685e-12  1.201376e-01  1.142421e+00 -1.398307e+00
              -6.735087e-03  1.467163e+00 -4.577402e-02  5.261200e-05", quiet = TRUE),
    nrow = 7)

  tau$A_1_5 <- matrix(
    scan(text = "5.952375e+00 -2.665165e+00  1.516589e+00  9.602737e-02 -6.547669e+00  1.721936e+00
              -1.555676e-01 -1.223167e+02  8.910711e+01 -2.137777e+01  5.352486e-01  1.634417e+02
              -4.421619e+01  3.702966e+00  1.077095e+02 -7.798154e+01  1.864825e+01 -4.615110e-01
              -1.437032e+02  3.918804e+01 -3.321332e+00  8.863768e+00 -6.471347e+00  1.565636e+00
              -4.037416e-02 -1.171173e+01  3.092645e+00 -2.503750e-01 -5.389123e+00  3.906715e+00
              -9.317209e-01  2.285102e-02  7.219457e+00 -1.979946e+00  1.692304e-01  1.632543e-01
              -1.235816e-01  3.163102e-02 -9.644010e-04 -2.089571e-01  4.926064e-02 -3.380261e-03
               3.616655e-05 -2.772220e-05  7.248774e-06 -2.365960e-07 -4.573620e-05  1.024572e-05
              -6.332864e-07 -3.489324e-02  2.518046e-02 -5.959302e-03  1.419892e-04  4.698732e-02
              -1.311024e-02  1.148966e-03 -7.681089e-06  5.534578e-06 -1.306550e-06  3.083095e-08
               1.035973e-05 -2.906102e-06  2.566389e-07 -1.673251e+02  1.226883e+02 -2.909396e+01
               7.343648e-01  2.249073e+02 -6.119645e+01  5.164918e+00", quiet = TRUE),
    nrow = 7)


  tau$A_5_10 <- rbind(
    c(-139.6445381,   425.6721706, -326.3024677, -27.36153163,  12.01878098, -3.655450e-02, -3.788487e-08,  1.554911e-03,  1.203641e-09,   562.9121367),
    c(87.6685893,  -243.6126993,  185.9780155,  15.80266620,  -6.82210340,  2.194062e-02,  2.284949e-08, -8.744362e-04, -6.761948e-10,  -320.1554409),
    c(-16.2933470,    42.7565168,  -32.3942750,  -2.80513070,   1.18311648, -4.092893e-03, -4.289413e-09,  1.496684e-04,  1.155655e-10,    56.1606779),
    c(0.3861494,    -0.4769648,  0.3529142,   0.03237058,  -0.01270732,  5.388348e-05,  5.741632e-11, -1.538782e-06, -1.182102e-12,    -0.6079357),
    c(251.1286638,  -856.5921761,  663.8547639,  54.31391421, -24.56950911,  6.727534e-02,  6.911057e-08, -3.229253e-03, -2.504112e-09, -1138.5802653),
    c(-184.1521518,   800.3905102, -632.4884072, -49.15903610,  23.66394490, -5.125947e-02, -5.157194e-08,  3.204510e-03,  2.492907e-09,  1077.7174769),
    c(183.6358765, -1033.1835539,  837.0245460,  60.69984342, -31.73480041,  4.938530e-02, 4.810099e-08, -4.449796e-03, -3.474150e-09, -1414.5947830)
  )

  tau$A_10_25 <- rbind(
    c(-334.0595182,  4.864037e+00,  -3.668889517, -0.2423163239,  0.1518617875, -5.561202e-04,  2.134014e-10,  2.278158e-05,  1.789987e-11,  5.188400e+00),
    c(192.7126831, -2.367611e+00,   2.013364346,  0.1277494909, -0.0769212844,  2.866288e-04, -1.089874e-10, -1.151519e-05, -9.050954e-12, -1.566204e+00),
    c(-33.5710167,  3.675872e-01,  -0.312245449, -0.0199548478,  0.0120090967, -4.445676e-05,  1.678260e-11,  1.793621e-06,  1.409718e-12,  4.577684e-01),
    c(0.5551671, -3.115692e-03,   0.002642717,  0.0001702175, -0.0001022234,  3.743850e-07, -1.390776e-13, -1.523990e-08, -1.197711e-14,  5.329921e-03),
    c(706.4932717, -1.178565e+01,   9.964299716,  0.6638424104, -0.3868855321,  1.464366e-03, -5.711352e-10, -5.799664e-05, -4.559005e-11, -1.162292e+01),
    c(-822.1302357,  1.790093e+01, -15.487163195, -0.9355838686,  0.6069634570, -2.630154e-03,  1.055902e-09,  9.254206e-05,  7.285793e-11,  1.985570e+01),
    c(2570.5029208, -1.648411e+02, 139.681907453,  8.9583503460, -5.4247388562,  1.264336e-02, -8.421040e-09, -8.087057e-04, -6.353406e-10, -2.152581e+02)
  )

  ztau <- qnorm(tau$TAU)
  tau1 <- 1 - tau$tau
  log_tau <- log(tau$tau)
  log_1_tau <- log(tau1)
  tauI <- 1/tau$tau
  tauI1 <- 1/tau1

  tau$B <- cbind(1, -log_tau, -log(tau1), log_tau^2, log(tau1)^2, -tauI, tauI^2, -tauI1, tauI1^2, ztau[1,])
  tau$z <- ztau

  tau
}

# Giles' algorithm 955
qpoisC.955 <- function(z, lambda){
  rlambda <- sqrt(lambda)
  q <- lambda + rlambda*z + (1 + 0.5*z^2)/3 - 1/rlambda/36*(z + 0.5*z^3)
  q - 0.0218/(q + 0.065*lambda)
}

# Frumento's approximation
qpoisC.me <- function(log.lambda, A, B){
  lambda1 <- exp(-log.lambda)
  U <- cbind(1, log.lambda, log.lambda^2, log.lambda^4, lambda1, lambda1^2, lambda1^4)
  temp <- (U %*% A) %*% t(B)
  if(any(temp < 0)) {temp[temp < 0] <- 0}
  temp
}

# Continuous Poisson distribution function: at integer y, ppoisC(y, lambda) = ppois(y - 1, lambda).
ppoisC <- function(y, lambda) pgamma(lambda, shape = y, scale = 1, lower.tail = FALSE)

# Very slow bisection (for lambda < 0.1)
qpoisC.bisec <- function(tau, lambda){

  yL <- qpois(tau, lambda)
  yR <- yL + 1

  for(i in 1:30){
    yC <- (yL + yR)/2
    tauC <- ppoisC(yC, lambda)
    wL <- (tauC < tau)
    wR <- !wL
    yL[wL] <- yC[wL]
    yR[wR] <- yC[wR]
  }
  (yL + yR)/2
}

# This function is for internal use only. The argument "obj" is created by Qpois.
qpoisC <- function(obj){

  lambda <- obj$lambda
  log.lambda <- obj$log.lambda
  w1 <- obj$w1
  w2 <- obj$w2
  w3 <- obj$w3
  w4 <- obj$w4
  w5 <- obj$w5
  w6 <- obj$w6
  tau <- obj$tau

  # I use different approximations depending on lambda.

  Q1 <- matrix(NA, length(lambda), ncol(tau$TAU))
  if(any(w1)){Q1[w1,] <- qpoisC.bisec(tau$TAU[w1,, drop = FALSE], lambda[w1])}
  if(any(w2)){Q1[w2,] <- qpoisC.me(log.lambda[w2], tau$A_0.1_1, tau$B)}
  if(any(w3)){Q1[w3,] <- qpoisC.me(log.lambda[w3], tau$A_1_5, tau$B)}
  if(any(w4)){Q1[w4,] <- qpoisC.me(log.lambda[w4], tau$A_5_10, tau$B)}
  if(any(w5)){Q1[w5,] <- qpoisC.me(log.lambda[w5], tau$A_10_25, tau$B)}
  if(any(w6)){Q1[w6,] <- qpoisC.955(tau$z[w6,, drop = FALSE], lambda[w6])}
  Q1
}

#######################################################################################################
# Estimating equation. Q(tau | x) depends on theta only through eta = x'theta, so that
# Q_theta = Q_eta x and Q_theta_theta = Q_eta_eta x x'. Q_eta and Q_eta_eta are computed by
# central differences in eta (with the same approximation of the quantile function), and
# F(y | x) by inverting the interpolant of Q on the grid, as in the generic engine of Qest.
# (In Qest <= 1.0.2, F(y | x) was computed exactly while A(F(y | x)) was interpolated on the
# approximated quantiles, which made the estimating equation non-smooth.)
QestPois.ee.u <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X")
  n <- nrow(X)
  tg <- tau$tau
  h <- 1e-4

  q0 <- Q(theta, tau, X)
  qshift <- function(delta){
    obj <- q0
    obj$log.lambda <- q0$log.lambda + delta
    obj$lambda <- exp(obj$log.lambda)
    qpoisC(obj)
  }
  Q0 <- q0$Q
  Qp <- qshift(h); Qm <- qshift(-h)
  Qe <- (Qp - Qm)/(2*h)
  wt <- rep(tau$wtau, each = n)
  B1 <- intblocks(Qe*wt, tg)
  ly <- locate(y, Q0, tg)
  g.i <- c(B1$Ab1 - interp(B1$A, ly))*X
  g <- c(w %*% g.i)
  fy <- ifelse(ly$out, 0, 1/ly$slope)
  Fy <- list(Fy = ly$F, delta = y - Q0)

  if(J){
    Qee <- (Qp - 2*Q0 + Qm)/h^2
    B2 <- intblocks(Qee*wt, tg)
    Fe <- ifelse(ly$out, 0, -interp(Qe, ly)/ly$slope)
    dg <- B2$Ab1 - interp(B2$A, ly) - islope(B1$A, ly)*Fe
    J <- crossprod(X*(w*dg), X)
  }

  list(g = g, g.i = g.i, J = J, Q1 = q0, Fy = Fy, fy = fy)
}

#######
# Quantile function

#' @rdname Qfamily
#' @export
Qpois <- function(offset = NULL) {
  Q <- function(theta, tau, X){

    offset <- attr(X, "offset")
    log.lambda <- c(X %*% theta + offset)
    lambda <- poisson()$linkinv(log.lambda)

    lambdac <- cut(lambda, c(0,.1,1,5,10,25,Inf), include.lowest = TRUE)
    levc <- levels(lambdac)
    w1 <- which(lambdac == levc[1])
    w2 <- which(lambdac == levc[2])
    w3 <- which(lambdac == levc[3])
    w4 <- which(lambdac == levc[4])
    w5 <- which(lambdac == levc[5])
    w6 <- which(lambdac == levc[6])

    obj <- list(lambda = lambda, log.lambda = log.lambda, w1 = w1, w2 = w2, w3 = w3,
                w4 = w4, w5 = w5, w6 = w6, #w7 = w7,
                tau = tau)
    obj$Q <- qpoisC(obj)
    obj

    # list(Q = Q1, lambda = lambda, tau = tau)
  }
  scale.pois <- function(X, y, z) {
    nms <- colnames(X)
    if(!("(Intercept)" %in% nms)) stop("You must include an intercept")
    Xc <- scale(X[,-1], center = TRUE, scale = TRUE)
    mX <- attributes(Xc)[[3]]; sX <- attributes(Xc)[[4]]
    Xc <- if(ncol(X) == 1) X else cbind("(Intercept)" = X[,1], Xc)
    sy <- 1; yc <- y /sy; zc <- z /sy
    Stats <- list(X = X, y = y, z = z, mX = mX, sX = sX, sy = sy)
    list(Xc = Xc, yc = yc, zc = zc, Stats = Stats)
  }
  descale.pois <- function(theta, Stats) {
    npar <- length(theta)
    if(npar > 1) {
      theta[2:npar] <- theta[2:npar] / Stats$sX
      theta[1] <- theta[1] - sum(theta[2:npar] * Stats$mX)
    }
    theta
  }
  initialize <- function(x, z, y, d, w, Q, start, tau, data, ok, Stats){
    if(any(y == 0)) warning("Please consider using jittered response (i.e., y + runif(n)). \nSee 'Qfamily' documentation for more details.")
    if(missing(start)) {
      y <- floor(y)
      use <- y < quantile(y, probs = 0.9)
      temp <- glm.fit(x[use, ], y[use], w[use], offset=attr(x, "offset")[use], family=poisson())
      theta <- temp$coefficients
    }
    else {
      npar <- length(ok)
      if(length(start) != npar) stop("Wrong size of 'start'")
      # if(any(is.na(start))) stop("NAs are not allowed in 'start'")
      start[is.na(start)] <- 0.
      theta <- start[ok]
      if(npar > 1){
        theta[2:npar] <- theta[2:npar]*Stats$sX
        theta[1] <- theta[1] + sum(theta[2:npar]*Stats$mX/Stats$sX)
      }
      nms <- colnames(x)
      names(theta) <- nms
    }
    theta
  }
  structure(list(family = poisson(), Q = Q, scale = scale.pois, descale = descale.pois,
                 offset = offset, initialize = initialize), class = "Qfamily")
}

#######################################################################################################
