# Gamma family: shape = exp(theta[1]), scale = exp(x'theta[-1]).

# Given an evaluation of a(tau) over a grid "tau" with half steps "dtau",
# returns A(tau) = int_0^tau a on the grid (Atau), at p (Ap), and at 1 (A1).
# The integral on (tau_K, 1) is approximated by a rectangle (in Qest <= 1.0.2 it
# was ignored, which biased the estimates by about 1-2% of their standard error).
# If int = FALSE, it simply returns a(p).
findAp.gamma <- function(atau, tau, dtau, p, int = TRUE){

  Atau <- Ap <- NULL
  n <- nrow(atau)
  A1 <- double(ncol(atau))
  for(j in 1:ncol(atau)){
    ap <- atau[,j]
    if(int){
      a <- c(0,ap)
      A <- cumsum((a[1:n] + a[2:(n + 1)])*dtau)
      A1[j] <- A[n] + (1 - tau[n])*ap[n]
      Ap <- cbind(Ap, approx(c(tau, 1), c(A, A1[j]), xout = p, rule = 2)$y)
    }
    else{
      A <- ap
      Ap <- cbind(Ap, approx(tau, A, xout = p, rule = 2)$y)
    }
    Atau <- cbind(Atau, A)
  }
  list(Atau = Atau, Ap = Ap, A1 = A1)
}

#######################################################################################################

#######################################################################################################
# First derivatives of Q(tau) w.r.t. the parameters.
# Both dQ.a and dQ.b are vectors of ntau elements.

derQtheta.gamma <- function(Q){

  a <- Q$a; b <- Q$b
  y <- Q$Q; tau <- Q$tau

  # a
  eps <- 1e-5
  y1 <- qgamma(tau, shape = a*exp(-eps), scale = 1)
  dQ.a <- (y - y1)/eps # To be multiplied by b[i]

  # b
  fy <- dgamma(y, shape = a, scale = 1)
  G1 <- a*(pgamma(y, shape = a + 1, scale = 1) - tau)
  dQ.b <- -G1/fy # To be multiplied by b[i]*X[i,j]


  cbind(dQ.a = dQ.a, dQ.b = dQ.b)
}

#######################################################################################################
# Second derivatives

der2Qtheta.gamma <- function(Q, Qtheta){

  a <- Q$a; b <- Q$b
  y <- Q$Q; tau <- Q$tau
  dQ.a <- Qtheta[,"dQ.a"]
  dQ.b <- Qtheta[,"dQ.b"]
  eps <- 1e-5

  # (a,a)
  y1 <- qgamma(tau, shape = a*exp(eps), scale = 1)
  dQ.a1 <- (y1 - y)/(eps)
  dQ.aa <- (dQ.a1 - dQ.a)/eps # To be multiplied by b[i]

  # (b,b)
  dQ.bb <- dQ.b # To be multiplied by b[i]*X[i,j]*X[i,h]

  # (a,b)
  dQ.ab <- dQ.a # To be multiplied by b[i]*X[i,j]

  cbind(dQ.aa = dQ.aa, dQ.ab = dQ.ab, dQ.bb = dQ.bb)
}

#######################################################################################################
# Estimating equation
QestGamma.ee.u <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X") # oppure X = data? Vedi tu
  n <- nrow(X)
  ntau <- tau$ntau

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, tau$tau, X) # Q(tau | x)/b, a vector of ntau elements
    Qtheta <- derQtheta.gamma(Q1) # 2 vectors of ntau elements
    Fy <- pgamma(y/Q1$b, shape = Q1$a, scale = 1)

    ap <- Qtheta*tau$wtau # a_theta(p) = wtau(p)*Qtheta(p)
    Ap <- findAp.gamma(ap, tau$tau, tau$dtau, Fy) # A_theta(p)
    Ay <- Ap$Ap # A_theta(F(y | x)) # matrix(n*2)
    AA1 <- findAp.gamma(Ap$Atau, tau$tau, tau$dtau, 1)$Ap # AA_theta(1) - two elements

    Ay <- Ay*Q1$b
    Ay <- cbind(Ay[,1], Ay[,2]*X)
    AA1 <- Q1$b%*%AA1
    AA1 <- cbind(AA1[,1], AA1[,2]*X)

    g.i <- AA1 - Ay
    g <- c(w %*% g.i)
    fy <- NULL
  }
  else{g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1; Qtheta <- EE$Qtheta; Fy <- EE$Fy; fy <- EE$fy}

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy))
    fy <- dgamma(y, shape = Q1$a, scale = Q1$b)
    XX <- tensorX(X)

    # A_theta_theta(p)
    Qthetatheta <- der2Qtheta.gamma(Q1, Qtheta)
    bp <- Qthetatheta*tau$wtau # a_theta_theta(p) = wtau(p)*Qthetatheta(p)
    Bp <- findAp.gamma(bp, tau$tau, tau$dtau, Fy) # A_theta_theta(p)
    By <- Bp$Ap # A_theta_theta(F(y)) # matrix(n,3)
    BB1 <- findAp.gamma(Bp$Atau, tau$tau, tau$dtau, 1)$Ap # AA_theta_theta(1) # - 3 elements
    Qtheta_Fy <- findAp.gamma(Qtheta, tau$tau, tau$dtau, Fy, int = FALSE)$Ap # Qheta(F(y)) # matrix(n,2)

    By <- By*Q1$b
    By <- cbind(By[,1], By[,2]*X, By[,3]*XX)
    BB1 <- Q1$b%*%BB1
    BB1 <- cbind(BB1[,1], BB1[,2]*X, BB1[,3]*XX)
    Qtheta_Fy <- Qtheta_Fy*Q1$b
    Qtheta_Fy <- cbind(Qtheta_Fy[,1], Qtheta_Fy[,2]*X)

    # Assemble

    npar <- length(theta)

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

  list(g = g, g.i = g.i, J = J, Q1 = Q1, Qtheta = Qtheta, Fy = Fy, fy = fy)
}

QestGamma.ee.c <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X")
  n <- nrow(X)
  ntau <- tau$ntau

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, tau$tau, X) # Q(tau | x)/b, a vector of ntau elements
    Qtheta <- derQtheta.gamma(Q1) # 2 vectors of ntau elements
    Fy <- pgamma(y/Q1$b, shape = Q1$a, scale = 1) # F(y | x)
    Sy <- 1 - Fy # S(y | x)

    ap <- Qtheta*tau$wtau # a_theta(p) = wtau(p)*Qtheta(p)
    Ap <- findAp.gamma(ap, tau$tau, tau$dtau, Fy) # A_theta(p)
    Ay <- Ap$Ap # A_theta(F(y | x)) # matrix(n*2)

    AAp <- findAp.gamma(Ap$Atau, tau$tau, tau$dtau, Fy) # AA_theta(p)
    AAy <- AAp$Ap # AA_theta(F(y | x)) # matrix(n*2)
    AA1 <- rbind(AAp$A1) # AA_theta(1) # two elements

    ###################

    Ay <- Ay*Q1$b
    Ay <- cbind(Ay[,1], Ay[,2]*X)

    AAy <- AAy*Q1$b
    AAy <- cbind(AAy[,1], AAy[,2]*X)

    AA1 <- Q1$b%*%AA1
    AA1 <- cbind(AA1[,1], AA1[,2]*X)

    ###################

    g.i <- (AA1 - d*Ay) + (1 - d)/Sy*(AAy - AA1)
    g <- c(w %*% g.i)
    fy <- NULL
  }
  else{
    g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1; Qtheta <- EE$Qtheta
    Fy <- EE$Fy; Sy <- EE$Sy; Ay <- EE$Ay; AAy <- EE$AAy; AA1 <- EE$AA1
  }

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy))
    fy <- dgamma(y, shape = Q1$a, scale = Q1$b)
    XX <- tensorX(X)

    # A_theta_theta(p)
    Qthetatheta <- der2Qtheta.gamma(Q1, Qtheta)
    bp <- Qthetatheta*tau$wtau # a_theta_theta(p) = wtau(p)*Qthetatheta(p)
    Bp <- findAp.gamma(bp, tau$tau, tau$dtau, Fy) # A_theta_theta(p)
    By <- Bp$Ap # A_theta_theta(F(y)) # matrix(n,3)
    BBp <- findAp.gamma(Bp$Atau, tau$tau, tau$dtau, Fy) # AA_theta_theta(p)
    BBy <- BBp$Ap # AA_theta_theta(F(y | x)) # matrix(n*3)
    BB1 <- rbind(BBp$A1) # AA_theta_theta(1) # 3 elements
    Qtheta_Fy <- findAp.gamma(Qtheta, tau$tau, tau$dtau, Fy, int = FALSE)$Ap # Qheta(F(y | x)) # matrix(n*2)

    ###################

    By <- By*Q1$b
    By <- cbind(By[,1], By[,2]*X, By[,3]*XX)

    BBy <- BBy*Q1$b
    BBy <- cbind(BBy[,1], BBy[,2]*X, BBy[,3]*XX)

    BB1 <- Q1$b%*%BB1
    BB1 <- cbind(BB1[,1], BB1[,2]*X, BB1[,3]*XX)

    Qtheta_Fy <- Qtheta_Fy*Q1$b
    Qtheta_Fy <- cbind(Qtheta_Fy[,1], Qtheta_Fy[,2]*X)

    ###################

    # Assemble

    npar <- length(theta)
    Uy <- (1 - d)/Sy

    h0 <- (BB1 - d*By) + Uy*(BBy - BB1)
    h0 <- .colSums(w*h0, n, npar*(npar + 1)/2)

    H0 <- matrix(NA, npar, npar)
    H0[lower.tri(H0, diag = TRUE)] <- h0; H0 <- t(H0)
    H0[lower.tri(H0, diag = TRUE)] <- h0

    H1 <- -Qtheta_Fy*(w*fy)
    H2 <- -d*wy*Qtheta_Fy + Uy*(Ay + (AAy - AA1)/Sy)
    H1 <- crossprod(H1, H2)

    # Finish

    H <- H0 + H1
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Q1 = Q1, Qtheta = Qtheta,
       Fy = Fy, fy = fy, Sy = Sy, Ay = Ay, AAy = AAy, AA1 = AA1)
}

QestGamma.ee.ct <- function(theta, eps, z, y, d, Q, w, data, tau, J = FALSE, EE){

  X <- attr(Q, "X") # oppure X = data? Vedi tu
  n <- nrow(X)
  ntau <- tau$ntau

  # Gradient
  if(missing(EE)){

    Q1 <- Q(theta, tau$tau, X) # Q(tau | x)/b, a vector of ntau elements
    Qtheta <- derQtheta.gamma(Q1) # 2 vectors of ntau elements

    Fy <- pgamma(y/Q1$b, shape = Q1$a, scale = 1) # F(y | x)
    Sy <- 1 - Fy # S(y | x)
    Fz <- pgamma(z/Q1$b, shape = Q1$a, scale = 1) # F(z | x)
    Sz <- 1 - Fz # S(z | x)

    ap <- Qtheta*tau$wtau # a_theta(p) = wtau(p)*Qtheta(p)
    Ap <- findAp.gamma(ap, tau$tau, tau$dtau, Fy) # A_theta(p)
    Ay <- Ap$Ap # A_theta(F(y | x)) # matrix(n*2)
    Az <- findAp.gamma(Ap$Atau, tau$tau, tau$dtau, Fz, int = FALSE)$Ap # A_theta(F(z | x)) # matrix(n*2)

    AAp <- findAp.gamma(Ap$Atau, tau$tau, tau$dtau, Fy) # AA_theta(p)
    AAy <- AAp$Ap # AA_theta(F(y | x)) # matrix(n*2)
    AAz <- findAp.gamma(AAp$Atau, tau$tau, tau$dtau, Fz, int = FALSE)$Ap # AA_theta(F(z | x)) # matrix(n*2)
    AA1 <- rbind(AAp$A1) # AA_theta(1) - two elements

    ###################

    Ay <- Ay*Q1$b
    Ay <- cbind(Ay[,1], Ay[,2]*X)

    Az <- Az*Q1$b
    Az <- cbind(Az[,1], Az[,2]*X)

    AAy <- AAy*Q1$b
    AAy <- cbind(AAy[,1], AAy[,2]*X)

    AAz <- AAz*Q1$b
    AAz <- cbind(AAz[,1], AAz[,2]*X)

    AA1 <- Q1$b%*%AA1
    AA1 <- cbind(AA1[,1], AA1[,2]*X)

    ###################

    gy.i <- -d*Ay + (1 - d)/Sy*(AAy - AA1)
    gz.i <- -1/Sz*(AAz - AA1)
    g.i <- gy.i + gz.i
    g <- c(w %*% g.i)
    fy <- fz <- NULL
  }
  else{
    g <- EE$g; g.i <- EE$g.i; Q1 <- EE$Q1; Qtheta <- EE$Qtheta
    Fy <- EE$Fy; Sy <- EE$Sy; Ay <- EE$Ay; AAy <- EE$AAy
    Fz <- EE$Fz; Sz <- EE$Sz; Az <- EE$Az; AAz <- EE$AAz
    AA1 <- EE$AA1
  }

  # Hessian

  if(J){

    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy))
    fy <- dgamma(y, shape = Q1$a, scale = Q1$b)
    fz <- dgamma(z, shape = Q1$a, scale = Q1$b)
    XX <- tensorX(X)

    # A_theta_theta(p)
    Qthetatheta <- der2Qtheta.gamma(Q1, Qtheta)
    bp <- Qthetatheta*tau$wtau # a_theta_theta(p) = wtau(p)*Qthetatheta(p)
    Bp <- findAp.gamma(bp, tau$tau, tau$dtau, Fy) # A_theta_theta(p)
    By <- Bp$Ap # A_theta_theta(F(y)) # matrix(n*3)
    Bz <- findAp.gamma(Bp$Atau, tau$tau, tau$dtau, Fz, int = FALSE)$Ap # A_theta_theta(F(z | x)) # matrix(n*3)
    BBp <- findAp.gamma(Bp$Atau, tau$tau, tau$dtau, Fy) # AA_theta_theta(p)
    BBy <- BBp$Ap # AA_theta_theta(F(y | x)) # matrix(n*3)
    BBz <- findAp.gamma(BBp$Atau, tau$tau, tau$dtau, Fz, int = FALSE)$Ap # AA_theta_theta(F(z | x)) # matrix(n*3)
    BB1 <- rbind(BBp$A1) # AA_theta_theta(1) # 3 elements
    Qtheta_Fy <- findAp.gamma(Qtheta, tau$tau, tau$dtau, Fy, int = FALSE)$Ap # Qheta(F(y | x)) # matrix(n*2)
    Qtheta_Fz <- findAp.gamma(Qtheta, tau$tau, tau$dtau, Fz, int = FALSE)$Ap # Qheta(F(z | x)) # matrix(n*2)

    ###################

    By <- By*Q1$b
    By <- cbind(By[,1], By[,2]*X, By[,3]*XX)

    Bz <- Bz*Q1$b
    Bz <- cbind(Bz[,1], Bz[,2]*X, Bz[,3]*XX)

    BBy <- BBy*Q1$b
    BBy <- cbind(BBy[,1], BBy[,2]*X, BBy[,3]*XX)

    BBz <- BBz*Q1$b
    BBz <- cbind(BBz[,1], BBz[,2]*X, BBz[,3]*XX)

    BB1 <- Q1$b%*%BB1
    BB1 <- cbind(BB1[,1], BB1[,2]*X, BB1[,3]*XX)

    Qtheta_Fy <- Qtheta_Fy*Q1$b
    Qtheta_Fy <- cbind(Qtheta_Fy[,1], Qtheta_Fy[,2]*X)

    Qtheta_Fz <- Qtheta_Fz*Q1$b
    Qtheta_Fz <- cbind(Qtheta_Fz[,1], Qtheta_Fz[,2]*X)

    ###################

    # Assemble

    npar <- length(theta)
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

    H2y <- -d*wy*Qtheta_Fy + Uy*(Ay + (AAy - AA1)/Sy)
    H2z <- -Uz*(Az + (AAz - AA1)/Sz)

    H1y <- crossprod(H1y, H2y)
    H1z <- crossprod(H1z, H2z)


    # Finish

    H <- H0 + H1y + H1z
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the hessian
  }

  list(g = g, g.i = g.i, J = J, Q1 = Q1, Qtheta = Qtheta,
       Fy = Fy, Sy = Sy, fy = fy, Ay = Ay, AAy = AAy,
       Fz = Fz, Sz = Sz, fz = fz, Az = Az, AAz = AAz,
       AA1 = AA1)
}

#######################################################################################################
# Quantile function

# Parametrization: gamma(shape = a, scale = b), with
 # a = exp(theta[1]), b = exp(X*theta[-1]).
# All quantities are computed using the fact that:
 # Q(tau | a,b)/b is a constant;
 # F(y | a,b) = F(y/b | a, 1).

# Q returns quantiles on the "tau" grid for a gamma(a,1).
# To get the actual quantiles of the i-th observation,
# the $Q value must be multiplied by b[i].

# This function does not receive a n*ntau grid of quantiles,
# but simply a vector of ntau quantiles.

#' @rdname Qfamily
#' @export
Qgamma <- function() {
  Q <- function(theta, tau, X){

    log.a <- theta[1]
    log.b <- c(X %*% theta[-1])
    a <- exp(log.a) # shape
    b <- Gamma("log")$linkinv(log.b) # scale

    Q1 <- qgamma(tau, shape = a, scale = 1)
    list(Q = Q1, a = a, b = b, tau = tau)
  }
  scale.gamma <- function(X, y, z) {
    nms <- colnames(X)
    if(!("(Intercept)" %in% nms)) stop("You must include an intercept")
    Xc <- scale(X[,-1], center = TRUE, scale = TRUE)
    mX <- attributes(Xc)[[3]]; sX <- attributes(Xc)[[4]]
    Xc <- if(ncol(X) == 1) X else cbind("(Intercept)" = X[,1], Xc)
    sy <- sd(y); yc <- y /sy; zc <- z /sy
    Stats <- list(X = X, y = y, z = z, mX = mX, sX = sX, sy = sy)
    list(Xc = Xc, yc = yc, zc = zc, Stats = Stats)
  }
  descale.gamma <- function(theta, Stats) {
    npar <- length(theta)
    if(npar == 2) {
      theta[2] <- theta[2] + log(Stats$sy)
    }
    else {
      theta[3:npar] <- theta[3:npar] / Stats$sX
      theta[2] <- theta[2] - sum(theta[3:npar] * Stats$mX) + log(Stats$sy)
    }
    theta
  }
  initialize <- function(x, z, y, d, w, Q, start, tau, data, ok, Stats){
    if(missing(start)) {
      use <- y < quantile(y, probs = 0.9)
      temp <- glm.fit(x[use, ], y[use], w[use], offset=attr(x, "offset")[use], family=Gamma("log"))
      temp$dispersion <- sum((w[use] * temp$residuals^2)[w[use] > 0]) / temp$df.residual
      names(temp$dispersion) <- "log(shape)"
      names(temp$coefficients) <- colnames(x)
      temp$coefficients[1] <- temp$coefficients[1] + log(temp$dispersion)
      theta <- c(-log(temp$dispersion), temp$coefficients)
    }
    else{
      npar <- length(ok) + 1
      if(length(start) != npar) stop("Wrong size of 'start'")
      # if(any(is.na(start))) stop("NAs are not allowed in 'start'")
      start[is.na(start)] <- 0.
      theta <- double(sum(ok) + 1)
      theta[1] <- start[1]
      theta[-1] <- start[-1][ok]
      if(npar == 2){theta[2] <- theta[2] - log(Stats$sy)}
      else{
        npar <- sum(ok) + 1
        theta[3:npar] <- theta[3:npar]*Stats$sX
        theta[2] <- theta[2] + sum(theta[3:npar]*Stats$mX/Stats$sX) - log(Stats$sy)
      }
      nms <- c("log(shape)", colnames(x))
      names(theta) <- nms
    }
    theta
  }
  fix.Fy.gamma <- function(fit, theta, tau, y, X, Q){ # Chiaramente, theta, y e X non standardizzati
    Q1 <- Q(theta, tau$tau, X)
    b <- Q1$b; Q1 <- Q1$Q
    Q1 <- t(matrix(Q1, length(Q1), length(b)))*b
    list(Fy = fit$Fy, delta = y - Q1)
  }
  structure(list(family = Gamma("log"), Q = Q, scale = scale.gamma, descale = descale.gamma,
                 fix.Fy = fix.Fy.gamma, initialize = initialize), class = "Qfamily")
}


#######################################################################################################
