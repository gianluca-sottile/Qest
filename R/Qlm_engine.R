# Building blocks of Qlm (Sottile and Frumento, 2022, Appendix D1):
# Q(tau | x) = x'beta + sigma * z(tau), z = qnorm, theta = (log(sigma), beta).
# W(tau) = int_0^tau w, ZW(tau) = int_0^tau w z, and their integrals on (0, 1)
# are computed once on a fine grid (4999 points in [1e-6, 1 - 1e-6]).

Qlm.bfun <- function(wtau, ...){

  r <- 4999
  p <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = r),2,2)
  dp <- p[-1] - p[-r]

  z <- qnorm(p)
  w <- wtau(p, ...)
  zw <- z*w
  W <- num.fun(dp,w)
  ZW <- num.fun(dp,zw)
  Wbar <- num.fun(dp,W)[r]
  ZWbar <- num.fun(dp,ZW)[r]

  z <- splinefun(p,z, method = "hyman") # faster than qnorm, no trouble at tau = 0 and tau = 1
  list(
    z = z,
    w = wtau,
    W = splinefun(p,W, method = "hyman"),
    ZW = approxfun(p,ZW, rule = 2), # may not be monotone, use approxfun
    Wbar = Wbar, ZWbar = ZWbar,
    opt = list(...)
  )
}

# Transform the scale of y and X in a linear model.
scalevars.Qlm <- function(X,y){

  sX <- apply(X, 2, sd)
  int <- which(sX == 0)
  vars <- which(sX > 0)
  is.int <- (length(int) > 0)
  mX <- colMeans(X)*is.int # zero if the model has no intercept

  X[,vars] <- scale(X[,vars], center = mX[vars], scale = sX[vars])

  my <- min(y); My <- max(y)
  y <- (y - my*is.int)/(My - my)*10 # y goes in (0,10). If no intercept, y is only scaled, not shifted
  Stats <- list(mX = mX, sX = sX, int = int, is.int = is.int, vars = vars, my = my, My = My)

  list(X = X, y = y, Stats = Stats)
}

# Transform back the coefficients of a linear model.
descalecoef.Qlm <- function(theta, Stats){

  beta <- theta[-1]
  log.sigma <- theta[1]

  beta[Stats$int] <- beta[Stats$int] - sum((beta*Stats$mX/Stats$sX)[Stats$vars])
  beta[Stats$vars] <- (beta/Stats$sX)[Stats$vars]

  beta <- beta * (Stats$My - Stats$my) / 10
  beta[Stats$int] <- beta[Stats$int] + Stats$my*Stats$is.int

  log.sigma <- log.sigma + log((Stats$My - Stats$my)/10)
  c(log.sigma, beta)
}

# Loss, gradient and Hessian.
qlmLoss <- function(theta, y, X, w, bfun){

  # building blocks

  sigma <- exp(theta[1]); beta <- theta[-1]
  delta <- c(y - X%*%beta)
  p.i <- pmax(pmin(pnorm(delta,0,sigma), 1-1e-6), 1e-6)
  deltaW <- bfun$W(p.i) - bfun$Wbar
  deltaZW <- bfun$ZW(p.i) - bfun$ZWbar

  Loss.i <- delta*(deltaW) - sigma*(deltaZW)
  Loss <- sum(w*Loss.i)
  return(list(Loss = Loss, p.i = p.i, delta = delta, deltaW = deltaW, deltaZW = deltaZW))
}

Qlm.ee.u <- function(theta, X, w, bfun, EE, J = FALSE){

  sigma <- exp(theta[1]); beta <- theta[-1]
  p.i <- EE$p.i; delta <- EE$delta; deltaW <- EE$deltaW; deltaZW <- EE$deltaZW

  # gradient

  g.sigma.i <- -sigma*deltaZW
  g.sigma <- sum(w*g.sigma.i)
  g.beta.i <- -(X*deltaW)
  g.beta <- c(w%*%g.beta.i)
  g <- c(g.sigma, g.beta)
  g.i <- cbind(g.sigma.i, g.beta.i)

  H <- NULL
  if(J){
    # hessian
    f <- dnorm(delta,0,sigma)
    wf <- w*f
    zi <- bfun$z(p.i)
    wi <- do.call(bfun$w, Ltau(bfun$opt, p.i))
    zwi <- zi*wi

    a <- wf*X*wi
    Hbeta <- t(a)%*%X
    Hsigma <- sigma^2*sum(wf*zwi*zi) + g[1]
    Hbetasigma <- sigma*rep(1,nrow(X))%*%(a*zi)
    H <- cbind(c(Hsigma, Hbetasigma), rbind(Hbetasigma, Hbeta))
  }

  # output

  list(g.i = g.i, g = g, H = H)
}
