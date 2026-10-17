# Normal family: Q(tau | x) = x'beta + sigma * qnorm(tau), theta = (log(sigma), beta).
# This does not actually return the quantile function, which is never needed during estimation.
#' @rdname Qfamily
#' @export
Qnorm <- function() {
  Q <- function(theta, X){
    sigma <- exp(theta[1])
    mu <- c(X %*% theta[-1])
    list(sigma = sigma, mu = mu, X = X)
  }
  scale.norm <- function(X, y, z) {
    sX <- apply(X, 2, sd)
    int <- which(sX == 0)
    vars <- which(sX > 0)
    is.int <- (length(int) > 0)
    mX <- colMeans(X)*is.int # zero if the model has no intercept

    Xc <- X; Xc[,vars] <- scale(X[,vars], center = mX[vars], scale = sX[vars])

    my <- min(y); My <- max(y)
    yc <- (y - my*is.int)/(My - my)*10 # y goes in (0,10). If no intercept, y is only scaled, not shifted
    zc <- (z - my*is.int)/(My - my)*10 # z goes in (0,10). If no intercept, y is only scaled, not shifted

    Stats <- list(X = X, y = y, z = z, mX = mX, sX = sX,
                  int = int, is.int = is.int, vars = vars, my = my, My = My)
    list(Xc = Xc, yc = yc, zc = zc, Stats = Stats)
  }
  descale.norm <- function(theta, Stats) {
    log.sigma <- theta[1]; beta <- theta[-1]
    K <- (Stats$My - Stats$my) / 10
    log.sigma <- log.sigma + log(K)

    beta[Stats$int] <- beta[Stats$int] - sum((beta * Stats$mX / Stats$sX)[Stats$vars])
    beta[Stats$vars] <- (beta / Stats$sX)[Stats$vars]

    beta <- beta * K
    beta[Stats$int] <- beta[Stats$int] + Stats$my * Stats$is.int

    c(log.sigma, beta)
  }
  bfun.norm <- function(wtau = NULL, wtauoptions = list()){

    if(is.null(wtau)){wtau <- function(tau, ...){rep(1,length(tau))}}

    r <- 4999
    p <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = r),2,2)
    dp <- p[-1] - p[-r]

    z <- qnorm(p)
    w <- callwtau(wtau, p, wtauoptions)
    zw <- z*w
    W <- num.fun(dp,w)
    ZW <- num.fun(dp,zw)
    iW <- num.fun(dp,W)
    iZW <- num.fun(dp,ZW)
    Wbar <- iW[r]
    ZWbar <- iZW[r]

    # I use hyman for monotone functions, approxfun otherwise.
    z <- splinefun(p,z, method = "hyman") # faster than qnorm, no trouble at tau = 0 and tau = 1
    list(
      z = z,
      w = wtau,
      W = splinefun(p,W, method = "hyman"),
      iW = splinefun(p,iW, method = "hyman"),
      ZW = approxfun(p,ZW, rule = 2),
      iZW = approxfun(p,iZW, rule = 2),
      Wbar = Wbar, ZWbar = ZWbar
    )
  }
  fix.Fy.norm <- function(fit, theta, tau, y, X, Q){
    Q1 <- Q(theta, X)
    Q1 <- qnorm(tau$TAU, Q1$mu, Q1$sigma)
    list(Fy = fit$Fy, delta = y - Q1)
  }
  initialize <- function(x, z, y, d, w, Q, start, tau, data, ok, Stats){
    if(missing(start)) {
      p1 <- .7; p2 <- .3
      temp <- rq.fit.br2(x, y, tau = .5)
      temp$dispersion <- (quantile(temp$residuals, p1) - quantile(temp$residuals, p2))/(qnorm(p1) - qnorm(p2))
      names(temp$dispersion) <- "log(sigma)"
      names(temp$coefficients) <- colnames(x)
      theta <- c(log(temp$dispersion), temp$coefficients)
    }
    else{
      npar <- length(ok) + 1
      if(length(start) != npar) stop("Wrong size of 'start'")
      # if(any(is.na(start))) stop("NAs are not allowed in 'start'")
      start[is.na(start)] <- 0.

      K <- (Stats$My - Stats$my) / 10
      theta <- double(sum(ok) + 1)
      theta[1] <- start[1] - log(K) #ok
      theta[-1] <- start[-1][ok]

      if(length(Stats$vars) == 0){theta[2] <- (theta[2] - Stats$my) / K}
      else{
        beta <- theta[-1]
        beta[Stats$vars] <- (beta * Stats$sX)[Stats$vars] / K
        if(Stats$is.int) beta[Stats$int] <- (beta[Stats$int]  - Stats$my) / K + sum((beta * Stats$mX / Stats$sX)[Stats$vars])
        theta[-1] <- beta
      }

      nms <- c("log(sigma)", colnames(x))
      names(theta) <- nms
    }
    theta
  }
  structure(list(family = gaussian(), Q = Q, scale = scale.norm, descale = descale.norm,
                 bfun = bfun.norm, fix.Fy = fix.Fy.norm, initialize = initialize),
            class = "Qfamily")
}


QestNorm.ee.u <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X")
  bfun <- tau$bfun
  n <- length(y)

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, X) # this is not Q(tau | x), actually
    Fy <- pnorm(y, Q1$mu, Q1$sigma)
    Fy <- pmin(pmax(1e-6, Fy), 1 - 1e-6)

    Ay.sigma <- Q1$sigma*bfun$ZW(Fy)
    Ay.beta <- bfun$W(Fy)
    AA1.sigma <- Q1$sigma*bfun$ZWbar
    AA1.beta <- bfun$Wbar

    g.i <- cbind(AA1.sigma - Ay.sigma, (AA1.beta - Ay.beta)*X)
    g <- c(w %*% g.i)
    fy <- NULL
  }
  else{g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1; Fy <- EE$Fy;
    Ay.sigma <- EE$Ay.sigma; AA1.sigma <- EE$AA1.sigma}

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy)) # w(F(y))
    fy <- dnorm(y, Q1$mu, Q1$sigma) # f(y)

    # Derivatives of Q w.r.t. theta, evaluated at F(y)
    Qtheta_Fy <- cbind(Q1$sigma*bfun$z(Fy), X)

    # A_theta_theta
    npar <- length(theta)
    By <- BB1 <- matrix(0, n, npar*(npar + 1)/2)
    By[,1] <- Ay.sigma
    BB1[,1] <- AA1.sigma

    # Assemble
    h0 <- BB1 - By
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)
    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    H1 <- -Qtheta_Fy*(w*fy)
    H2 <- -wy*Qtheta_Fy
    H1 <- crossprod(H1, H2)

    # Finish

    H <- H0 + H1
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Q1 = Q1, Fy = Fy, fy = fy, Ay.sigma = Ay.sigma, AA1.sigma = AA1.sigma)
}

QestNorm.ee.c <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X")
  bfun <- tau$bfun
  n <- length(y)

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, X) # this is not Q(tau | x), actually
    Fy <- pnorm(y, Q1$mu, Q1$sigma)
    Fy <- pmin(pmax(1e-6, Fy), 1 - 1e-6)
    Sy <- 1 - Fy

    Ay.sigma <- Q1$sigma*bfun$ZW(Fy)
    Ay.beta <- bfun$W(Fy)
    AAy.sigma <- Q1$sigma*bfun$iZW(Fy)
    AAy.beta <- bfun$iW(Fy)
    AA1.sigma <- Q1$sigma*bfun$ZWbar
    AA1.beta <- bfun$Wbar

    g.i.sigma <- (AA1.sigma - d*Ay.sigma) + (1 - d)/Sy*(AAy.sigma - AA1.sigma)
    g.i.beta <- (AA1.beta - d*Ay.beta) + (1 - d)/Sy*(AAy.beta - AA1.beta)
    g.i <- cbind(g.i.sigma, g.i.beta*X)

    g <- c(w %*% g.i)
    fy <- NULL
  }
  else{
    g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1; Fy <- EE$Fy; Sy <- EE$Sy
    Ay.sigma <- EE$Ay.sigma; AAy.sigma <- EE$AAy.sigma; AA1.sigma <- EE$AA1.sigma
    Ay.beta <- EE$Ay.beta; AAy.beta <- EE$AAy.beta; AA1.beta <- EE$AA1.beta
  }

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy)) # w(F(y))
    fy <- dnorm(y, Q1$mu, Q1$sigma) # f(y)

    # Derivatives of Q w.r.t. theta, evaluated at F(y)
    Qtheta_Fy <- cbind(Q1$sigma*bfun$z(Fy), X)

    # A_theta_theta
    npar <- length(theta)
    By <- BBy <- BB1 <- matrix(0, n, npar*(npar + 1)/2)
    By[,1] <- Ay.sigma
    BBy[,1] <- AAy.sigma
    BB1[,1] <- AA1.sigma

    # Assemble

    Uy <- (1 - d)/Sy
    h0 <- (BB1 - d*By) + Uy*(BBy - BB1)
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)

    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    Ay <- cbind(Ay.sigma, Ay.beta*X)
    AAy_AA1 <- cbind(AAy.sigma - AA1.sigma, (AAy.beta - AA1.beta)*X)
    H1 <- -Qtheta_Fy*(w*fy)
    H2 <- -d*wy*Qtheta_Fy + Uy*(Ay + (AAy_AA1)/Sy)
    H1 <- crossprod(H1, H2)

    # Finish

    H <- H0 + H1
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(
    g = g, g.i = g.i, J = J, Q1 = Q1, Fy = Fy, fy = fy, Sy = Sy,
    Ay.sigma = Ay.sigma, AAy.sigma = AAy.sigma, AA1.sigma = AA1.sigma,
    Ay.beta = Ay.beta, AAy.beta = AAy.beta, AA1.beta = AA1.beta
  )
}

QestNorm.ee.ct <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X")
  bfun <- tau$bfun
  n <- length(y)

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, X) # this is not Q(tau | x), actually

    Fy <- pnorm(y, Q1$mu, Q1$sigma)
    Fy <- pmin(pmax(1e-6, Fy), 1 - 1e-6)
    Sy <- 1 - Fy

    Fz <- pnorm(z, Q1$mu, Q1$sigma)
    Fz <- pmin(pmax(1e-6, Fz), 1 - 1e-6)
    Sz <- 1 - Fz

    Ay.sigma <- Q1$sigma*bfun$ZW(Fy)
    Ay.beta <- bfun$W(Fy)
    Az.sigma <- Q1$sigma*bfun$ZW(Fz)
    Az.beta <- bfun$W(Fz)
    AAy.sigma <- Q1$sigma*bfun$iZW(Fy)
    AAy.beta <- bfun$iW(Fy)
    AAz.sigma <- Q1$sigma*bfun$iZW(Fz)
    AAz.beta <- bfun$iW(Fz)
    AA1.sigma <- Q1$sigma*bfun$ZWbar
    AA1.beta <- bfun$Wbar


    gy.i.sigma <- -d*Ay.sigma + (1 - d)/Sy*(AAy.sigma - AA1.sigma)
    gz.i.sigma <- -1/Sz*(AAz.sigma - AA1.sigma)
    gy.i.beta <- -d*Ay.beta + (1 - d)/Sy*(AAy.beta - AA1.beta)
    gz.i.beta <- -1/Sz*(AAz.beta - AA1.beta)

    g.i <- cbind(gy.i.sigma + gz.i.sigma, (gy.i.beta + gz.i.beta)*X)
    g <- c(w %*% g.i)
    fy <- fz <- NULL
  }
  else{
    g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1
    Fy <- EE$Fy; Sy <- EE$Sy
    Ay.sigma <- EE$Ay.sigma; AAy.sigma <- EE$AAy.sigma
    Ay.beta <- EE$Ay.beta; AAy.beta <- EE$AAy.beta
    Fz <- EE$Fz; Sz <- EE$Sz
    Az.sigma <- EE$Az.sigma; AAz.sigma <- EE$AAz.sigma
    Az.beta <- EE$Az.beta; AAz.beta <- EE$AAz.beta
    AA1.sigma <- EE$AA1.sigma; AA1.beta <- EE$AA1.beta
  }

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy)) # w(F(y))
    fy <- dnorm(y, Q1$mu, Q1$sigma) # f(y)
    fz <- dnorm(z, Q1$mu, Q1$sigma) # f(z)

    # Derivatives of Q w.r.t. theta, evaluated at F(y) and F(z)
    Qtheta_Fy <- cbind(Q1$sigma*bfun$z(Fy), X)
    Qtheta_Fz <- cbind(Q1$sigma*bfun$z(Fz), X)

    # A_theta_theta
    npar <- length(theta)

    By <- BBy <- Bz <- BBz <- BB1 <- matrix(0, n, npar*(npar + 1)/2)
    By[,1] <- Ay.sigma
    BBy[,1] <- AAy.sigma
    Bz[,1] <- Az.sigma
    BBz[,1] <- AAz.sigma
    BB1[,1] <- AA1.sigma

    # Assemble

    Uy <- (1 - d)/Sy
    Uz <- 1/Sz

    h0y <- -d*By + Uy*(BBy - BB1)
    h0z <- -Uz*(BBz - BB1)
    h0 <- h0z + h0y
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)

    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    H1y <- -Qtheta_Fy*(w*fy)
    H1z <- -Qtheta_Fz*(w*fz)


    Ay <- cbind(Ay.sigma, Ay.beta*X)
    AAy_AA1 <- cbind(AAy.sigma - AA1.sigma, (AAy.beta - AA1.beta)*X)
    Az <- cbind(Az.sigma, Az.beta*X)
    AAz_AA1 <- cbind(AAz.sigma - AA1.sigma, (AAz.beta - AA1.beta)*X)

    H2y <- -d*wy*Qtheta_Fy + Uy*(Ay + (AAy_AA1)/Sy)
    H2z <- -Uz*(Az + (AAz_AA1)/Sz)

    H1y <- crossprod(H1y, H2y)
    H1z <- crossprod(H1z, H2z)

    # Finish

    H <- H0 + H1y + H1z
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(
    g = g, g.i = g.i, J = J, Q1 = Q1,
    Fy = Fy, fy = fy, Sy = Sy,
    Ay.sigma = Ay.sigma, AAy.sigma = AAy.sigma,
    Ay.beta = Ay.beta, AAy.beta = AAy.beta,
    Fz = Fz, fz = fz, Sz = Sz,
    Az.sigma = Az.sigma, AAz.sigma = AAz.sigma,
    Az.beta = Az.beta, AAz.beta = AAz.beta,
    AA1.sigma = AA1.sigma, AA1.beta = AA1.beta
  )
}

