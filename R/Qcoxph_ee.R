# Qcoxph estimating equations and their Jacobian, with the closed forms of
# Sottile and Frumento (2022, Appendix D2). Building blocks below.
# The argument "eps" is not used, but is required for compatibility with Qest.newton.

# Estimating equation + Jacobian of a Cox model with censored data.
QCox.ee.c <- function(theta, eps, z, y, d, X, w, knots, tau, J = FALSE, EE){

  n <- nrow(X)
  q <- ncol(X)
  k <- length(knots) - 2 # n. of internal knots
  Ltau <- function(opt, tau){opt$tau <- tau; opt}

  #############################################################################
  #############################################################################
  # Estimating equation #######################################################
  #############################################################################
  #############################################################################

  if(missing(EE)){

    # Building blocks
    BBy <- coxBB(theta, y, X, knots, tau)

    # A_beta(p) - to be A*X
    Ab <- A_beta_fun(BBy)
    A_beta <- Ab$A
    AA_beta <- Ab$AA
    AA1_beta <- Ab$AA1

    # A_gamma(p)
    Ag <- A_gamma_fun(BBy)
    A_gamma <- Ag$A
    AA_gamma <- Ag$AA
    AA1_gamma <- Ag$AA1

    #############################################################################
    # Estimating equation #######################################################
    #############################################################################

    Ay <- cbind(A_beta*X, A_gamma) # A_theta
    AAy <- cbind(AA_beta*X, AA_gamma) # AA_theta
    AA1 <- cbind(AA1_beta*X, AA1_gamma) # AA1_theta
    g.i <- -d*Ay + (1 - d)/BBy$Sy*(AAy - AA1) + AA1

    # Regularization. If F(y) = 1, the estimating equation has a 0/0 = NaN.

    outLy <- BBy$outL; outRy <- BBy$outR
    g.i[outLy,] <- (d*AA1)[outLy,]
    g.i[outRy,] <- (AA1 - d*Ay)[outRy,]
    g <- c(w %*% g.i)
  }
  else{
    g <- EE$g; g.i <- EE$g.i
    BBy <- EE$BBy; Ay <- EE$Ay; AAy <- EE$AAy; AA1 <- EE$AA1
  }

  #############################################################################
  #############################################################################
  # Jacobian ##################################################################
  #############################################################################
  #############################################################################

  if(J){

    # A_beta_beta(p) - to be t(A*X)%*%X
    Abb <- A_beta_beta_fun(BBy)

    # A_gamma_jh(p)
    if(k > 0){Agmix <- A_gamma_gamma_mix_fun(BBy)}

    # A_gamma_jj(p)
    Agg <- A_gamma_gamma_fun(BBy)

    # A_beta_gamma(p) - to be t(x)%*%A to get a matrix q*(k + 1)
    Abg <- A_beta_gamma_fun(BBy)


    #############################################################################
    #############################################################################
    # Assemble the Jacobian #####################################################
    #############################################################################
    #############################################################################

    genH <- function(A, d, U, outL, outR){
      H <- -d*A$A + U*(A$AA - A$AA1) + A$AA1
      H <- cbind(H)
      H[outLy,] <- d[outLy]*cbind(A$AA1)[outLy,]
      H[outRy,] <- cbind(A$AA1)[outRy,] - d[outRy]*cbind(A$A)[outRy,]
      H
    }

    Xw <- X*w
    Sy <- BBy$Sy; Fy <- BBy$Fy; fy <- BBy$fy
    U <- (1 - d)/Sy
    outLy <- BBy$outL; outRy <- BBy$outR

    H_beta <- t(c(genH(Abb, d, U, outLy, outRy))*X)%*%Xw
    H_gamma_mix <- (if(k > 0) .colSums(genH(Agmix, d, U, outLy, outRy)*w, n, k*(k + 1)/2) else NULL)
    H_gamma_jj <- .colSums(genH(Agg, d, U, outLy, outRy)*w, n, k + 1)
    H_beta_gamma <- t(Xw)%*%genH(Abg, d, U, outLy, outRy)


    # Assemble

    H0 <- matrix(NA, q + k + 1, q + k + 1)
    i_beta <- (1:q); i_gamma <- (q + 1):(q + k + 1)
    H0[i_beta, i_beta] <- H_beta
    H0[i_beta, i_gamma] <- H_beta_gamma; H0[i_gamma, i_beta] <- t(H_beta_gamma)
    H0gamma <- matrix(0, k + 1, k + 1)
    if(k > 0){H0gamma[lower.tri(H0gamma)] <- H_gamma_mix; H0gamma <- H0gamma + t(H0gamma)}
    diag(H0gamma) <- H_gamma_jj
    H0[i_gamma, i_gamma] <- H0gamma

    fy[c(outLy, outRy)] <- 0
    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy))
    H1 <- -fy*BBy$Qtheta
    H2 <- -d*wy*BBy$Qtheta + U*(Ay + (AAy - AA1)/Sy)

    # Finish

    H <- crossprod(w*H1, H2) + H0
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the jacobian
  }

  list(g = g, g.i = g.i, J = J,
       Fy = BBy$Fy, fy = BBy$fy,
       Fz = NULL, fz = NULL,
       BBy = BBy, Ay = Ay, AAy = AAy,
       BBz = NULL, Az = NULL, AAz = NULL,
       AA1 = AA1
  )
}

# Estimating equation + Jacobian of a Cox model with censored and truncated data.
QCox.ee.ct <- function(theta, eps, z, y, d, X, w, knots, tau, J = FALSE, EE){

  n <- nrow(X)
  q <- ncol(X)
  k <- length(knots) - 2 # n. of internal knots
  Ltau <- function(opt, tau){opt$tau <- tau; opt}

  #############################################################################
  #############################################################################
  # Estimating equation #######################################################
  #############################################################################
  #############################################################################

  if(missing(EE)){

    # Building blocks
    BBy <- coxBB(theta, y, X, knots, tau)
    BBz <- coxBB(theta, z, X, knots, tau)

    # A_beta(p) - to be A*X
    Aby <- A_beta_fun(BBy); Abz <- A_beta_fun(BBz)
    A_betay <- Aby$A; A_betaz <- Abz$A
    AA_betay <- Aby$AA; AA_betaz <- Abz$AA
    AA1_beta <- Aby$AA1

    # A_gamma(p)
    Agy <- A_gamma_fun(BBy); Agz <- A_gamma_fun(BBz)
    A_gammay <- Agy$A; A_gammaz <- Agz$A
    AA_gammay <- Agy$AA; AA_gammaz <- Agz$AA
    AA1_gamma <- Agy$AA1

    #############################################################################
    # Estimating equation #######################################################
    #############################################################################

    Sy <- BBy$Sy; Sz <- BBz$Sy; Fy <- BBy$Fy; Fz <- BBz$Fy

    Ay <- cbind(A_betay*X, A_gammay) # A_theta(Fy)
    AAy <- cbind(AA_betay*X, AA_gammay) # AA_theta(Fy)
    Az <- cbind(A_betaz*X, A_gammaz) # A_theta(Fz)
    AAz <- cbind(AA_betaz*X, AA_gammaz) # AA_theta(Fz)
    AA1 <- cbind(AA1_beta*X, AA1_gamma) # AA1_theta
    gy.i <- -d*Ay + (1 - d)/Sy*(AAy - AA1)
    gz.i <- -1/Sz*(AAz - AA1)

    # Regularization. If F(y) = 1 or F(z) = 1, the estimating equation has a 0/0 = NaN.

    outLy <- BBy$outL; outRy <- BBy$outR
    outLz <- BBz$outL; outRz <- BBz$outR
    gy.i[outLy,] <- (-(1 - d)*AA1)[outLy,]
    gz.i[outLz,] <- AA1[outLz,]
    gy.i[outRy,] <- (-d*Ay)[outRy,]
    gz.i[outRz,] <- 0

    g.i <- gy.i + gz.i
    g <- c(w %*% g.i)
  }
  else{
    g <- EE$g; g.i <- EE$g.i
    BBy <- EE$BBy; BBz <- EE$BBz
    Ay <- EE$Ay; AAy <- EE$AAy; Az <- EE$Az; AAz <- EE$AAz; AA1 <- EE$AA1
  }

  #############################################################################
  #############################################################################
  # Jacobian ##################################################################
  #############################################################################
  #############################################################################

  if(J){

    # A_beta_beta(p) - to be t(A*X)%*%X
    Abby <- A_beta_beta_fun(BBy)
    Abbz <- A_beta_beta_fun(BBz)

    # A_gamma_jh(p)
    if(k > 0){
      Agmixy <- A_gamma_gamma_mix_fun(BBy)
      Agmixz <- A_gamma_gamma_mix_fun(BBz)
    }

    # A_gamma_jj(p)
    Aggy <- A_gamma_gamma_fun(BBy)
    Aggz <- A_gamma_gamma_fun(BBz)

    # A_beta_gamma(p) - to be t(x)%*%A to get a matrix q*(k + 1)
    Abgy <- A_beta_gamma_fun(BBy)
    Abgz <- A_beta_gamma_fun(BBz)

    #############################################################################
    #############################################################################
    # Assemble the Jacobian #####################################################
    #############################################################################
    #############################################################################

    genH <- function(Ay, Az, d, Uy, Uz, outLy, outLz, outRy, outRz){

      Hy <- -d*Ay$A + Uy*(Ay$AA - Ay$AA1)
      Hz <- -Uz*(Az$AA - Az$AA1)
      Hy <- cbind(Hy); Hz <- cbind(Hz)

      Hy[outLy,] <- -(1 - d[outLy])*cbind(Ay$AA1)[outLy,]
      Hz[outLz,] <- cbind(Az$AA1)[outLz,]
      Hy[outRy,] <- -d[outRy]*cbind(Ay$A)[outRy,]
      Hz[outRz,] <- 0

      Hy + Hz
    }

    Xw <- X*w
    Sy <- BBy$Sy; Sz <- BBz$Sy; Fy <- BBy$Fy; Fz <- BBz$Fy; fy <- BBy$fy; fz <- BBz$fy
    Uy <- (1 - d)/Sy; Uz <- 1/Sz
    outLy <- BBy$outL; outRy <- BBy$outR
    outLz <- BBz$outL; outRz <- BBz$outR

    H_beta <- t(c(genH(Abby, Abbz, d, Uy, Uz, outLy, outLz, outRy, outRz))*X)%*%Xw
    H_gamma_mix <- (if(k > 0) .colSums(genH(Agmixy, Agmixz, d, Uy, Uz, outLy, outLz, outRy, outRz)*w, n, k*(k + 1)/2) else NULL)
    H_gamma_jj <- .colSums(genH(Aggy, Aggz, d, Uy, Uz, outLy, outLz, outRy, outRz)*w, n, k + 1)
    H_beta_gamma <- t(Xw)%*%genH(Abgy, Abgz, d, Uy, Uz, outLy, outLz, outRy, outRz)


    # Assemble

    H0 <- matrix(NA, q + k + 1, q + k + 1)
    i_beta <- (1:q); i_gamma <- (q + 1):(q + k + 1)
    H0[i_beta, i_beta] <- H_beta
    H0[i_beta, i_gamma] <- H_beta_gamma; H0[i_gamma, i_beta] <- t(H_beta_gamma)
    H0gamma <- matrix(0, k + 1, k + 1)
    if(k > 0){H0gamma[lower.tri(H0gamma)] <- H_gamma_mix; H0gamma <- H0gamma + t(H0gamma)}
    diag(H0gamma) <- H_gamma_jj
    H0[i_gamma, i_gamma] <- H0gamma

    fy[c(outLy, outRy)] <- 0
    fz[c(outLz, outRz)] <- 0
    wy <- do.call(tau$wfun, Ltau(tau$opt, Fy))
    H1y <- -fy*BBy$Qtheta
    H2y <- -d*wy*BBy$Qtheta + Uy*(Ay + (AAy - AA1)/Sy)
    H1z <- -fz*BBz$Qtheta
    H2z <- Uz*(Az + (AAz - AA1)/Sz)

    # Finish

    H <- crossprod(w*H1y, H2y) - crossprod(w*H1z, H2z) + H0
    J <- H # So if J = FALSE, return J = FALSE; and if J = TRUE, return the jacobian
  }

  list(g = g, g.i = g.i, J = J,
       Fy = BBy$Fy, fy = BBy$fy,
       Fz = BBz$Fy, fz = BBz$fy,
       BBy = BBy, Ay = Ay, AAy = AAy,
       BBz = BBz, Az = Az, AAz = AAz,
       AA1 = AA1
  )
}



# Building blocks for Cox model

coxBB <- function(theta, y, X, knots, tau){

  n <- nrow(X)
  q <- ncol(X)
  k <- length(knots) - 2 # n. of internal knots
  beta <- (theta[1:q]) # Cox regression coefficients
  egamma <- exp(theta[(q + 1):length(theta)]) # coefficients of H0

  delta <- function(obj, e, k){t(t(obj[,2:(k + 2)] - obj[,1:(k + 1)])/e)}
  deltaj <- function(obj, e, b, k){o <- list(); for(j in 1:(k + 1)){o[[j]] <- delta(t(t(obj)*(e[j]*b[,j])), e, k)}; o}

  #############################################################################
  #############################################################################
  #############################################################################

  eta <- c(X%*%beta); HR <- exp(eta); HR1 <- 1/HR

  # Building blocks - part 1
  by <- plfcox(y, knots, deriv = 1) # b(y)
  H0y <- c(by$b%*%egamma) # H0(y)
  Fystar <- eta + log(H0y) # log(-log(1 - F(y | x)))
  outL <- which(Fystar <= tau$taustarL); outR <- which(Fystar >= tau$taustarR)
  Fystar <- pmax(pmin(Fystar, tau$taustarR), tau$taustarL)
  Hy <- exp(Fystar) # H(y | x)
  Sy <- exp(-Hy) # S(y | x)
  Fy <- 1 - Sy # F(y | x)


  # Building blocks - part 2
  bknots <- plfcox(knots, knots, deriv = 0) # b(knots)
  H0knots <- c(bknots%*%egamma) # H0(knots)
  Fknotstar <- eta + rep(1,n)%*%t(log(H0knots)) # log(-log(1 - F(knots | x)))
  Fknotstar[,k + 2] <- tau$taustarR # "as if" the last knot was Inf ("use all quantiles"). This is fundamental.
  Fknotstar <- pmax(pmin(Fknotstar, tau$taustarR), tau$taustarL)
  Hknots <- exp(Fknotstar) # H(knots | x)
  Sknots <- exp(-Hknots) # S(knots | x)
  Fknots <- 1 - Sknots # F(knots | x)


  # Building blocks - part 3 - Derivatives of the quantile function
  h0 <- c(by$b1 %*% egamma) # h0(y)
  hy <- h0*HR # h(y | x)
  fy <- hy*Sy # f(y | x)
  Qbeta <- -X*(Hy/hy)
  Qgamma <- -t(t(by$b/h0)*egamma)
  Qtheta <- cbind(Qbeta, Qgamma)


  # Building blocks - part 4

  ii <- (Fknotstar < Fystar)
  mFknotsp <- pmax0(Sknots - Sy)
  mFknots1 <- Sknots

  mins <- pmin(Fknots, Fy)
  minstar <- pmin(Fknotstar, Fystar)

  mins <- c(mins)
  minstar <- c(minstar)
  Fknots <- c(Fknots)
  Fknotstar <- c(Fknotstar)

  # w blocks

  w1mins <- matrix(tau$w1fun(mins), n, k + 2)
  wmins <- matrix(do.call(tau$wfun, Ltau(tau$opt, mins)), n, k + 2)
  Wmins <- matrix(tau$Wfun(mins), n, k + 2)
  iWmins <- matrix(tau$iWfun(mins), n, k + 2)

  w1knots <- matrix(tau$w1fun(Fknots), n, k + 2)
  w1knotsp <- w1knots*mFknotsp
  w1knots1 <- w1knots*mFknots1

  wknots <- matrix(do.call(tau$wfun, Ltau(tau$opt, Fknots)), n, k + 2)
  wknotsp <- wknots*mFknotsp
  wknots1 <- wknots*mFknots1

  Wknots <- matrix(tau$Wfun(Fknots), n, k + 2)
  Wknotsp <- Wknots*mFknotsp
  Wknots1 <- Wknots*mFknots1

  iWknots <- matrix(tau$iWfun(Fknots), n, k + 2)
  iWknotsp <- iWknots*mFknotsp

  wmins_2wknots <- wmins - 2*wknots
  Wmins_Wknots <- Wmins - Wknots
  w1knots1_wknots <- w1knots1 - wknots

  dWmins <- delta(Wmins, egamma, k)
  diWmins_dWknotsp <- delta(iWmins + Wknotsp, egamma, k)
  diWknots_dWknots1 <- delta(iWknots + Wknots1, egamma, k)


  # r blocks

  r1mins <- matrix(do.call(tau$r1fun, Ltau(tau$opt, mins)), n, k + 2)*HR1
  rmins <- matrix(do.call(tau$rfun, Ltau(tau$opt, minstar)), n, k + 2)*HR1
  Rmins <- matrix(tau$Rfun(minstar), n, k + 2)*HR1
  iRmins <- matrix(tau$iRfun(minstar), n, k + 2)*HR1

  r1knots <- matrix(do.call(tau$r1fun, Ltau(tau$opt, Fknots)), n, k + 2)*HR1
  r1knotsp <- r1knots*mFknotsp
  r1knots1 <- r1knots*mFknots1

  rknots <- matrix(do.call(tau$rfun, Ltau(tau$opt, Fknotstar)), n, k + 2)*HR1
  rknotsp <- rknots*mFknotsp
  rknots1 <- rknots*mFknots1

  Rknots <- matrix(tau$Rfun(Fknotstar), n, k + 2)*HR1
  Rknotsp <- Rknots*mFknotsp
  Rknots1 <- Rknots*mFknots1

  iRknots <- matrix(tau$iRfun(Fknotstar), n, k + 2)*HR1
  iRknotsp <- iRknots*mFknotsp

  rmins_2rknots <- rmins - 2*rknots
  Rmins_Rknots <- Rmins - Rknots
  r1knots1_rknots <- r1knots1 - rknots

  dRmins <- delta(Rmins, egamma, k)
  diRmins_dRknotsp <- delta(iRmins + Rknotsp, egamma, k)
  diRknots_dRknots1 <- delta(iRknots + Rknots1, egamma, k)


  # Building blocks with Vbeta

  vbeta <- t(t(Sknots)*c(H0knots))*HR
  Vbeta <- vbeta*ii

  dwminsVbeta <- delta(wmins*Vbeta, egamma, k)
  dwknots1vbeta <- delta(wknots1*vbeta, egamma, k)
  Owbeta <- delta(Wmins_Wknots*Vbeta + wknotsp*vbeta, egamma, k)

  drminsVbeta <- delta(rmins*Vbeta, egamma, k)
  drknots1vbeta <- delta(rknots1*vbeta, egamma, k)
  Orbeta <- delta(Rmins_Rknots*Vbeta + rknotsp*vbeta, egamma, k)


  # Building blocks with Vgamma
  # Note: in practice, I need Vgamma_hj = Vgamma*exp(gamma[j])*bj(knots)

  vgamma <- Sknots*HR
  Vgamma <- vgamma*ii

  dwminsVgamma <- deltaj(wmins*Vgamma, egamma, bknots, k)
  dwknots1vgamma <- deltaj(wknots1*vgamma, egamma, bknots, k)
  Owgamma <- deltaj(Wmins_Wknots*Vgamma + wknotsp*vgamma, egamma, bknots, k)

  drminsVgamma <- deltaj(rmins*Vgamma, egamma, bknots, k)
  drknots1vgamma <- deltaj(rknots1*vgamma, egamma, bknots, k)
  Orgamma <- deltaj(Rmins_Rknots*Vgamma + rknotsp*vgamma, egamma, bknots, k)


  #############################################################################
  # Output - only the strict necessary ########################################
  #############################################################################


  list(
   n = n, q = q, k = k,
   beta = beta, egamma = egamma, knots = knots,
   HR = HR, HR1 = HR1, fy = fy, Qtheta = Qtheta,
   by = by, Hy = Hy, Sy = Sy, Fy = Fy,
   bknots = bknots, H0knots = H0knots, Sknots = Sknots,
   outL = outL, outR = outR,
   ii = ii, c1 = knots[1:(k + 1)]*egamma - H0knots[1:(k + 1)],

   ###

   w1mins = w1mins, wmins = wmins, Wmins = Wmins,
   w1knotsp = w1knotsp, w1knots1 = w1knots1,
   wknots = wknots, wknotsp = wknotsp, wknots1 = wknots1,
   Wknots = Wknots,

   wmins_2wknots = wmins_2wknots,
   Wmins_Wknots = Wmins_Wknots,
   w1knots1_wknots = w1knots1_wknots,

   dWmins = dWmins,
   diWmins_dWknotsp = diWmins_dWknotsp,
   diWknots_dWknots1 = diWknots_dWknots1,

   ###

   r1mins = r1mins, rmins = rmins, Rmins = Rmins,
   r1knotsp = r1knotsp, r1knots1 = r1knots1,
   rknots = rknots, rknotsp = rknotsp, rknots1 = rknots1,
   Rknots = Rknots,

   rmins_2rknots = rmins_2rknots,
   Rmins_Rknots = Rmins_Rknots,
   r1knots1_rknots = r1knots1_rknots,

   dRmins = dRmins,
   diRmins_dRknotsp = diRmins_dRknotsp,
   diRknots_dRknots1 = diRknots_dRknots1,

   ###

   vbeta = vbeta,
   Vbeta = Vbeta,

   dwminsVbeta = dwminsVbeta,
   dwknots1vbeta = dwknots1vbeta,
   Owbeta = Owbeta,

   drminsVbeta = drminsVbeta,
   drknots1vbeta = drknots1vbeta,
   Orbeta = Orbeta,

   ###

   vgamma = vgamma,
   Vgamma = Vgamma,

   dwminsVgamma = dwminsVgamma,
   dwknots1vgamma = dwknots1vgamma,
   Owgamma = Owgamma,

   drminsVgamma = drminsVgamma,
   drknots1vgamma = drknots1vgamma,
   Orgamma = Orgamma
  )
}






# A_beta
A_beta_fun <- function(BB){

  n <- BB$n; k <- BB$k

  c1 <- BB$c1
  c2A <- BB$dRmins
  c2AA <- BB$diRmins_dRknotsp
  c2AA1 <- BB$diRknots_dRknots1

  #########################################

  A <- .rowSums(t(t(BB$dwminsVbeta)*c1) - c2A + BB$drminsVbeta, n, k + 1)
  AA <- .rowSums(t(t(BB$Owbeta)*c1) - c2AA + BB$Orbeta, n, k + 1)
  AA1 <- .rowSums(t(t(BB$dwknots1vbeta)*c1) - c2AA1 + BB$drknots1vbeta, n, k + 1)

  list(A = A, AA = AA, AA1 = AA1)
}




# A_gamma
A_gamma_fun <- function(BB){

  egamma <- BB$egamma
  n <- BB$n; k <- BB$k
  H0knots <- BB$H0knots[1:(k + 1)]

  # A ################################################################

  cAw <- BB$dWmins
  cAr <- BB$dRmins
  A <- t(t(cAw)*H0knots) - cAr

  # AA ################################################################

  cAAw <- BB$diWmins_dWknotsp
  cAAr <- BB$diRmins_dRknotsp
  AA <- t(t(cAAw)*H0knots) - cAAr

  # AA1 ################################################################

  cAA1w <- BB$diWknots_dWknots1
  cAA1r <- BB$diRknots_dRknots1
  AA1 <- t(t(cAA1w)*H0knots) - cAA1r

  if(k == 0){return(list(A = A, AA = AA, AA1 = AA1))}
  for(j in 1:k){ # When j = k + 1, the result is zero.

    ej <- egamma[j]*BB$bknots[,j]

    indh <- j:(k + 1) # Which "h" I will use for this "l"
    ejh <- ej[indh]
    egammah <- egamma[indh]
    c1h <- BB$c1[indh]

    # A ################################################################

    t1 <- t(t(BB$dwminsVgamma[[j]][,indh])*c1h)
    t2 <- (BB$drminsVgamma[[j]][,indh])
    t3 <- -t(t(BB$dWmins[,indh])*ejh)
    A[,j] <- A[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

    # AA ################################################################

    t1 <- t(t(BB$Owgamma[[j]][,indh])*c1h)
    t2 <- (BB$Orgamma[[j]][,indh])
    t3 <- -t(t(BB$diWmins_dWknotsp[,indh])*ejh)
    AA[,j] <- AA[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

    # AA1 ################################################################

    t1 <- t(t(BB$dwknots1vgamma[[j]][,indh])*c1h)
    t2 <- (BB$drknots1vgamma[[j]][,indh])
    t3 <- -t(t(BB$diWknots_dWknots1[,indh])*ejh)
    AA1[,j] <- AA1[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)
  }
  list(A = A, AA = AA, AA1 = AA1)
}






# A_beta_beta
A_beta_beta_fun <- function(BB){

  delta <- function(obj, egamma){k <- ncol(obj); t(t(obj[,2:k] - obj[,1:(k - 1)])/egamma)}
  n <- BB$n; k <- BB$k
  egamma <- BB$egamma
  c1 <- BB$c1

  # Preliminary quantities

  HH <- 1 - BB$HR%*%t(BB$H0knots)
  Vbeta <- BB$Vbeta; vbeta <- BB$vbeta
  Vbeta2 <- Vbeta^2
  vbeta2 <- vbeta^2
  Vbetastar <- Vbeta*HH
  vbetastar <- vbeta*HH

  # A ################################################################

  uw <- BB$w1mins*Vbeta2
  Uw <- BB$wmins*Vbetastar
  dw <- delta(uw + Uw, egamma)

  ur <- BB$r1mins*Vbeta2
  Ur <- BB$rmins*Vbetastar
  dr <- delta(ur + Ur, egamma)

  cr <- BB$dRmins - 2*BB$drminsVbeta

  A <- .rowSums(t(t(dw)*c1) + (dr + cr), n, k + 1)

  # AA ################################################################

  uw <- (BB$wmins_2wknots)*Vbeta2 + BB$w1knotsp*vbeta2
  Uw <- (BB$Wmins_Wknots)*Vbetastar + BB$wknotsp*vbetastar
  dw <- delta(uw + Uw, egamma)

  ur <- (BB$rmins_2rknots)*Vbeta2 + BB$r1knotsp*vbeta2
  Ur <- (BB$Rmins_Rknots)*Vbetastar + BB$rknotsp*vbetastar
  dr <- delta(ur + Ur, egamma)

  cr <- BB$diRmins_dRknotsp - 2*BB$Orbeta

  AA <- .rowSums(t(t(dw)*c1) + (dr + cr), n, k + 1)

  # AA1 ################################################################

  uw <- (BB$w1knots1_wknots)*vbeta2
  Uw <- BB$wknots1*vbetastar
  dw <- delta(uw + Uw, egamma)

  ur <- (BB$r1knots1_rknots)*vbeta2
  Ur <- BB$rknots1*vbetastar
  dr <- delta(ur + Ur, egamma)

  cr <- BB$diRknots_dRknots1 - 2*BB$drknots1vbeta

  AA1 <- .rowSums(t(t(dw)*c1) + (dr + cr), n, k + 1)

  #######################

  list(A = A, AA = AA, AA1 = AA1)
}





# A_(gamma_j, gamma_l) with l > j
A_gamma_gamma_mix_fun <- function(BB){

  delta <- function(obj, egamma){k <- ncol(obj); t(t(obj[,2:k] - obj[,1:(k - 1)])/egamma)}
  n <- BB$n; k <- BB$k
  egamma <- BB$egamma
  H0knots <- BB$H0knots[1:(k + 1)]

  A <- AA <- AA1 <- NULL
  for(j in 1:k){
    indl <- (j + 1):(k + 1) # Which "l" I will use for this "j"
    ej <- egamma[j]*BB$bknots[,j]

    # A ################################################################

    cjA.1 <- BB$dwminsVgamma[[j]]
    cjA.2 <- BB$dWmins
    cjA.3 <- BB$drminsVgamma[[j]]
    cjA <- t(t(cjA.1)*H0knots) + t(t(cjA.2)*ej[1:(k + 1)]) - cjA.3

    # AA ################################################################

    cjAA.1 <- BB$Owgamma[[j]]
    cjAA.2 <- BB$diWmins_dWknotsp
    cjAA.3 <- BB$Orgamma[[j]]
    cjAA <- t(t(cjAA.1)*H0knots) + t(t(cjAA.2)*ej[1:(k + 1)]) - cjAA.3

    # AA1 ################################################################

    cjAA1.1 <- BB$dwknots1vgamma[[j]]
    cjAA1.2 <- BB$diWknots_dWknots1
    cjAA1.3 <- BB$drknots1vgamma[[j]]
    cjAA1 <- t(t(cjAA1.1)*H0knots) + t(t(cjAA1.2)*ej[1:(k + 1)]) - cjAA1.3

    for(l in indl){

      A_jl <- cjA[,l]
      AA_jl <- cjAA[,l]
      AA1_jl <- cjAA1[,l]

      if(l < k + 1){ # When l = k + 1, the result is zero.

        indh <- l:(k + 1) # Which "h" I will use for this "l"
        indh2 <- l:(k + 2) # Again, including last
        el <- egamma[l]*BB$bknots[,l]

        v_hlj <- -BB$vgamma[,indh2]*(BB$HR%*%t(ej[indh2]*el[indh2]))
        v_hl_hj <- -v_hlj*BB$Sknots[,indh2]
        V_hlj <- v_hlj*BB$ii[,indh2]
        V_hl_hj <- v_hl_hj*BB$ii[,indh2]

        ejh <- ej[indh]; elh <- el[indh]
        egammah <- egamma[indh]
        c1h <- BB$c1[indh]

        # A ################################################################

        w_V_hl_hj <- BB$w1mins[,indh2]*V_hl_hj
        w_V_hlj <- BB$wmins[,indh2]*V_hlj

        t1 <- delta(w_V_hl_hj + w_V_hlj, egammah)
        t1 <- t(t(t1)*c1h)

        r_V_hl_hj <- BB$r1mins[,indh2]*V_hl_hj
        r_V_hlj <- BB$rmins[,indh2]*V_hlj
        t2 <- (r_V_hl_hj + r_V_hlj)
        t2 <- delta(t2, egammah)

        dw_mixj <- BB$dwminsVgamma[[j]][,indh]
        dw_mixl <- BB$dwminsVgamma[[l]][,indh]
        t3 <- -(t(t(dw_mixj)*elh) + t(t(dw_mixl)*ejh))

        A_jl <- A_jl + .rowSums(t1 + t2 + t3, n, k - l + 2)

        # AA ################################################################

        w_V_hl_hj <- (BB$wmins_2wknots)[,indh2]*V_hl_hj
        w_V_hlj <- (BB$Wmins_Wknots)[,indh2]*V_hlj
        w_v_hl_hj <- BB$w1knotsp[,indh2]*v_hl_hj
        w_v_hlj <- BB$wknotsp[,indh2]*v_hlj
        t1 <- delta(w_V_hl_hj + w_V_hlj + w_v_hl_hj + w_v_hlj, egammah)
        t1 <- t(t(t1)*c1h)

        r_V_hl_hj <- (BB$rmins_2rknots)[,indh2]*V_hl_hj
        r_V_hlj <- (BB$Rmins_Rknots)[,indh2]*V_hlj
        r_v_hl_hj <- BB$r1knotsp[,indh2]*v_hl_hj
        r_v_hlj <- BB$rknotsp[,indh2]*v_hlj
        t2 <- (r_V_hl_hj + r_V_hlj + r_v_hl_hj + r_v_hlj)
        t2 <- delta(t2, egammah)

        dw_mixj <- BB$Owgamma[[j]][,indh]
        dw_mixl <- BB$Owgamma[[l]][,indh]
        t3 <- -(t(t(dw_mixj)*elh) + t(t(dw_mixl)*ejh))

        AA_jl <- AA_jl + .rowSums(t1 + t2 + t3, n, k - l + 2)

        # AA1 ################################################################

        w_v_hl_hj <- (BB$w1knots1_wknots)[,indh2]*v_hl_hj
        w_v_hlj <- BB$wknots1[,indh2]*v_hlj
        t1 <- delta(w_v_hl_hj + w_v_hlj, egammah)
        t1 <- t(t(t1)*c1h)

        r_v_hl_hj <- (BB$r1knots1_rknots)[,indh2]*v_hl_hj
        r_v_hlj <- BB$rknots1[,indh2]*v_hlj
        t2 <- (r_v_hl_hj + r_v_hlj)
        t2 <- delta(t2, egammah)

        dw_mixj <- BB$dwknots1vgamma[[j]][,indh]
        dw_mixl <- BB$dwknots1vgamma[[l]][,indh]
        t3 <- -(t(t(dw_mixj)*elh) + t(t(dw_mixl)*ejh))

        AA1_jl <- AA1_jl + .rowSums(t1 + t2 + t3, n, k - l + 2)
      }

      A <- cbind(A, A_jl)
      AA <- cbind(AA, AA_jl)
      AA1 <- cbind(AA1, AA1_jl)
    }
  }

  list(A = A, AA = AA, AA1 = AA1)
}




# A_(gamma_j, gamma_j)
A_gamma_gamma_fun <- function(BB){

  delta <- function(obj, egamma){k <- ncol(obj); t(t(obj[,2:k] - obj[,1:(k - 1)])/egamma)}
  listj <- function(x){out <- NULL; for(j in 1:length(x)){out <- cbind(out, x[[j]][,j])}; out}
  n <- BB$n; k <- BB$k

  egamma <- BB$egamma
  H0knots <- BB$H0knots[1:(k + 1)]
  c1 <- BB$c1

  # Preliminary quantities

  e <- egamma*BB$bknots[k + 2,] # Final value reached by each b(y)

  vgammaHR <- BB$vgamma*BB$HR
  vgamma2 <- BB$vgamma^2
  VgammaHR <- BB$Vgamma*BB$HR
  Vgamma2 <- BB$Vgamma^2

  v <- t(t(BB$vgamma[,-1])*(e/egamma)) # V(j + 1, j)/exp(gamma)
  V <- t(t(BB$Vgamma[,-1])*(e/egamma)) # V(j + 1, j)/exp(gamma)

  # A ################################################################

  cAw <- 2*(V*BB$wmins[,-1]) - BB$dWmins
  cAr <- 2*(V*BB$rmins[,-1]) - BB$dRmins
  A <- t(t(cAw)*H0knots) - cAr

  # AA ################################################################

  cAAw <- 2*(listj(BB$Owgamma)) - BB$diWmins_dWknotsp
  cAAr <- 2*(listj(BB$Orgamma)) - BB$diRmins_dRknotsp
  AA <- t(t(cAAw)*H0knots) - cAAr

  # AA1 ################################################################

  cAA1w <- 2*(listj(BB$dwknots1vgamma)) - BB$diWknots_dWknots1
  cAA1r <- 2*(listj(BB$drknots1vgamma)) - BB$diRknots_dRknots1
  AA1 <- t(t(cAA1w)*H0knots) - cAA1r


  if(k == 0){return(list(A = A, AA = AA, AA1 = AA1))}
  for(j in 1:k){ # When j = k + 1, the result is zero.

    indh <- j:(k + 1) # Which "h" I will use for this "l"
    indh2 <- j:(k + 2) # Again, including last

    ej <- egamma[j]*BB$bknots[,j]
    v_hjj <- -t(t(vgammaHR[,indh2])*ej[indh2])
    v2_hj <- t(t(vgamma2[,indh2])*(ej[indh2]^2))
    V_hjj <- -t(t(VgammaHR[,indh2])*ej[indh2])
    V2_hj <- t(t(Vgamma2[,indh2])*(ej[indh2]^2))

    ejh <- ej[indh]
    egammah <- egamma[indh]
    c1h <- c1[indh]

    # A ################################################################

    w_V_hjj <- BB$wmins[,indh2]*V_hjj
    w_V2_hj <- BB$w1mins[,indh2]*V2_hj
    t1 <- delta(w_V_hjj + w_V2_hj, egammah)
    t1 <- t(t(t1)*c1h)

    r_V_hjj <- BB$rmins[,indh2]*V_hjj
    r_V2_hj <- BB$r1mins[,indh2]*V2_hj
    t2 <- (r_V_hjj + r_V2_hj)
    t2 <- delta(t2, egammah)

    dw_mixj <- 2*BB$dwminsVgamma[[j]][,indh] + BB$dWmins[,indh]
    t3 <- -t(t(dw_mixj)*ejh)

    A[,j] <- A[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

    # AA ################################################################

    w_V_hjj <- (BB$Wmins_Wknots[,indh2])*V_hjj
    w_V2_hj <- (BB$wmins_2wknots[,indh2])*V2_hj
    w_v_hjj <- BB$wknotsp[,indh2]*v_hjj
    w_v2_hj <- BB$w1knotsp[,indh2]*v2_hj
    t1 <- delta(w_V_hjj + w_V2_hj + w_v_hjj + w_v2_hj, egammah)
    t1 <- t(t(t1)*c1h)


    r_V_hjj <- (BB$Rmins_Rknots[,indh2])*V_hjj
    r_V2_hj <- (BB$rmins_2rknots[,indh2])*V2_hj
    r_v_hjj <- BB$rknotsp[,indh2]*v_hjj
    r_v2_hj <- BB$r1knotsp[,indh2]*v2_hj
    t2 <- (r_V_hjj + r_V2_hj + r_v_hjj + r_v2_hj)
    t2 <- delta(t2, egammah)

    dw_mixj <- 2*BB$Owgamma[[j]][,indh] + BB$diWmins_dWknotsp[,indh]
    t3 <- -t(t(dw_mixj)*ejh)

    AA[,j] <- AA[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

    # AA1 ################################################################

    w_v_hjj <- BB$wknots1[,indh2]*v_hjj
    w_v2_hj <- (BB$w1knots1_wknots[,indh2])*v2_hj
    t1 <- delta(w_v_hjj + w_v2_hj, egammah)
    t1 <- t(t(t1)*c1h)

    r_v_hjj <- BB$rknots1[,indh2]*v_hjj
    r_v2_hj <- (BB$r1knots1_rknots[,indh2])*v2_hj
    t2 <- (r_v_hjj + r_v2_hj)
    t2 <- delta(t2, egammah)

    dw_mixj <- 2*BB$dwknots1vgamma[[j]][,indh] + BB$diWknots_dWknots1[,indh]
    t3 <- -t(t(dw_mixj)*ejh)

    AA1[,j] <- AA1[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)
  }
  list(A = A, AA = AA, AA1 = AA1)
}







# A_beta_gamma
A_beta_gamma_fun <- function(BB){

  delta <- function(obj, egamma){k <- ncol(obj); t(t(obj[,2:k] - obj[,1:(k - 1)])/egamma)}
  n <- BB$n; q <- BB$q; k <- BB$k
  egamma <- BB$egamma
  H0knots <- BB$H0knots[1:(k + 1)]
  c1 <- BB$c1

  # A ################################################################

  cAw <- BB$dwminsVbeta
  cAr <- BB$dRmins - BB$drminsVbeta
  A <- t(t(cAw)*H0knots) + cAr

  # AA ################################################################

  cAAw <- BB$Owbeta
  cAAr <- BB$diRmins_dRknotsp - BB$Orbeta
  AA <- t(t(cAAw)*H0knots) + cAAr

  # AA1 ################################################################

  cAA1w <- BB$dwknots1vbeta
  cAA1r <- BB$diRknots_dRknots1 - BB$drknots1vbeta
  AA1 <- t(t(cAA1w)*H0knots) + cAA1r

  if(k == 0){return(list(A = A, AA = AA, AA1 = AA1))}
  for(j in 1:k){ # When j = k + 1, the result is zero.

    ej <- egamma[j]*BB$bknots[,j]

    indh <- j:(k + 1) # Which "h" I will use for this "l"
    indh2 <- j:(k + 2) # Again, including last

    v_mixj <- t(t(-BB$vbeta[,indh2]*BB$HR + BB$vgamma[,indh2])*(ej[indh2]))
    vv_j <- t(t(BB$vbeta[,indh2]*BB$vgamma[,indh2])*(ej[indh2]))
    V_mixj <- v_mixj*BB$ii[,indh2]
    VV_j <- vv_j*BB$ii[,indh2]

    ejh <- ej[indh]
    egammah <- egamma[indh]
    c1h <- c1[indh]

    # A ################################################################

    w_VV_j <- BB$w1mins[,indh2]*VV_j
    w_V_mixj <- BB$wmins[,indh2]*V_mixj
    t1 <- delta(w_VV_j + w_V_mixj, egammah)
    t1 <- t(t(t1)*c1h)

    r_VV_j <- BB$r1mins[,indh2]*VV_j
    r_V_mixj <- BB$rmins[,indh2]*V_mixj
    t2 <- (r_VV_j + r_V_mixj)
    t2 <- delta(t2, egammah)

    dw_j <- BB$dwminsVbeta[,indh]
    dr_j <- BB$drminsVgamma[[j]][,indh]
    t3 <- -t(t(dw_j)*ejh) - dr_j

    A[,j] <- A[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

    # AA ################################################################

    w_VV_j <- (BB$wmins[,indh2] - 2*BB$wknots[,indh2])*VV_j
    w_V_mixj <- (BB$Wmins[,indh2] - BB$Wknots[,indh2])*V_mixj
    w_vv_j <- BB$w1knotsp[,indh2]*vv_j
    w_v_mixj <- BB$wknotsp[,indh2]*v_mixj
    t1 <- delta(w_VV_j + w_V_mixj + w_vv_j + w_v_mixj, egammah)
    t1 <- t(t(t1)*c1h)

    r_VV_j <- (BB$rmins[,indh2] - 2*BB$rknots[,indh2])*VV_j
    r_V_mixj <- (BB$Rmins[,indh2] - BB$Rknots[,indh2])*V_mixj
    r_vv_j <- BB$r1knotsp[,indh2]*vv_j
    r_v_mixj <- BB$rknotsp[,indh2]*v_mixj
    t2 <- (r_VV_j + r_V_mixj + r_vv_j + r_v_mixj)
    t2 <- delta(t2, egammah)

    dw_j <- BB$Owbeta[,indh]
    dr_j <- BB$Orgamma[[j]][,indh]
    t3 <- -t(t(dw_j)*ejh) - dr_j

    AA[,j] <- AA[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

    # AA1 ################################################################

    w_vv_j <- (BB$w1knots1[,indh2] - BB$wknots[,indh2])*vv_j
    w_v_mixj <- BB$wknots1[,indh2]*v_mixj
    t1 <- delta(w_vv_j + w_v_mixj, egammah)
    t1 <- t(t(t1)*c1h)

    r_vv_j <- (BB$r1knots1[,indh2] - BB$rknots[,indh2])*vv_j
    r_v_mixj <- BB$rknots1[,indh2]*v_mixj
    t2 <- (r_vv_j + r_v_mixj)
    t2 <- delta(t2, egammah)

    dw_j <- BB$dwknots1vbeta[,indh]
    dr_j <- BB$drknots1vgamma[[j]][,indh]
    t3 <- -t(t(dw_j)*ejh) - dr_j

    AA1[,j] <- AA1[,j] + .rowSums(t1 + t2 + t3, n, k - j + 2)

  }
  list(A = A, AA = AA, AA1 = AA1)
}

# Basis of a piecewise linear function.
plfcox <- function(y, knots, deriv = 0){ # knots incl. boundary

  if(length(knots) == 2){
    if(deriv == 0){return(cbind(y - knots[1]))}
    else{return(list(b = cbind(y - knots[1]), b1 = cbind(rep(1,length(y)))))}
  }
  intknots <- knots[2:(length(knots) - 1)]
  k <- length(intknots)
  ind1 <- 1
  ind2 <- NULL
  for (j in intknots){ind1 <- cbind(ind1,(y > j))}
  for (j in k:1) {ind2 <- cbind(ind2, ind1[, j] - ind1[, j + 1])}
  ind2 <- cbind(ind1[, k + 1], ind2)[, (k + 1):1, drop = FALSE]
  ind1 <- cbind(ind1, 0)

  b <- NULL
  for (j in 1:(k + 1)){b <- cbind(b, (y - knots[j])*ind2[, j] + (knots[j + 1] - knots[j])*ind1[, j + 1])}
  colnames(b) <- paste("b", 1:(k + 1), sep = "")
  attr(b, "knots") <- knots
  if(deriv == 0){return(b)}
  else{return(list(b = b, b1 = ind2))}
}

# Transform the scale of y and X in a Cox model.
scalevars.Qcoxph <- function(X,z,y,knots){

  mX <- colMeans(X)
  sX <- apply(X, 2, sd)
  X <- scale(X, center = mX, scale = sX)

  My <- max(y)
  y <- y/My*10
  z <- z/My*10
  knots <- knots/My*10
  Stats <- list(mX = mX, sX = sX, My = My)

  list(X = X, z = z, y = y, knots = knots, Stats = Stats)
}

descalecoef.Qcoxph <- function(theta, Stats){

  q <- length(Stats$mX)
  beta <- theta[1:q]
  log.gamma <- theta[(q + 1):length(theta)]

  beta <- beta/Stats$sX
  log.gamma <- log.gamma + log(10/Stats$My) - sum(beta*Stats$mX)

  c(beta, log.gamma)
}

# Positive slopes for the starting points of H0.
adjust.coef <- function(theta){
  id <- which(theta <= 0)
  id1 <- which(theta > 0)
  k <- length(theta)
  if(length(id) != 0){
    for(i in rev(id)){
      theta[i] <- if(i == 1) theta[i+1]/k else median(theta[-unique(c(k:i, id))])/k
    }
  }
  theta
}

# Loss of Cox model.
# If the data are not censored and not truncated, there is an analytical expression:
# A <- .rowSums(t(t(BBy$dWmins)*BBy$c1) + BBy$dRmins, n, k + 1)
# AA1 <- .rowSums(t(t(BBy$diWknots_dWknots1)*BBy$c1) + BBy$diRknots_dRknots1, n, k + 1)
# L <- sum(w*(y*(tau$Wfun(Fy) - tau$iWfun(1)) + AA1 - A))
# However, I assume that cox model will (almost) always be used with censored or truncated data
# (which is why I did not bother to implement QCox.ee.u), and I keep the code simple by
# directly computing the numeric integral that is needed with censored or truncated data.
coxLoss <- function(theta, z, y, d, X, w, knots, tau, type, Fy, Fz){

  Ltau <- function(opt, tau){opt$tau <- tau; opt}
  n <- nrow(X)
  q <- ncol(X)
  k <- length(knots) - 2
  beta <- (theta[1:q]) # Cox regression coefficients
  egamma <- exp(theta[(q + 1):length(theta)]) # coefficients of H0
  HR1 <- exp(-c(X%*%beta))

  # Quantile function evaluated on a grid

  taugrid <- c(1e-6, (1:99)/100, 1 - 1e-6); dtau <- 0.01; ntau <- length(taugrid)
  TAU <- t(matrix(taugrid, ntau, n))
  log1_tau <- log(1 - taugrid)
  u <- -HR1*rep(log1_tau, each = n) # vector of length n*ntau

  knotsinv <- c(plfcox(knots, knots = knots)%*%egamma) # knots of the inverse of b
  knotsinv_long <- unique(c(0, knotsinv, -max(HR1)*log1_tau[ntau])) # plfcox can project outside range, approx cannot.
  binv <- plfcox(knotsinv_long, knotsinv)

  Qcoxtemp <- c(binv%*%(1/egamma))
  Qcox <- approx(knotsinv_long, Qcoxtemp, xout = u, method = "linear")$y + knots[1]
  qq <- matrix(Qcox, n, ntau) # A matrix n*ntau.

  # Pseudo-Loss for censored or truncated data. With "u" data, omega = omegay.

  delta <- y - qq
  omegay <- (delta <= 0)
  omega <- omegay + (1 - d)*omegay*(TAU - 1)/(1 - pmin(Fy, 1 - 1e-10))
  if(type == "ct"){
    deltaz <- z - qq
    omegaz <- (deltaz <= 0)
    omega <- omega - omegaz - omegaz*(TAU - 1)/(1 - pmin(Fz, 1 - 1e-10)) + TAU
  }

  wtau <- do.call(tau$wfun, Ltau(tau$opt, taugrid))
  out <- delta*(TAU - omega)
  out <- out %*% (dtau*wtau)
  sum(w*out)
}
