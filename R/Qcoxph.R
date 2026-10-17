#' Q-Estimation of Proportional Hazards Regression Models
#'
#' Fits a proportional hazards model with the Q-estimator of Sottile and
#' Frumento (2022). The baseline cumulative hazard \eqn{H_0(t)} is piecewise
#' linear, so that the model is fully parametric, and
#' \eqn{Q(\tau \mid x) = H_0^{-1}(-\exp(-x^\top\beta)\log(1 - \tau))}.
#'
#' @param formula a formula with a `survival::Surv()` response: `Surv(time, event)`
#'   for right-censored data, or `Surv(start, stop, event)` for left-truncated
#'   and right-censored data. Times must be positive.
#' @param weights an optional vector of weights.
#' @param start optional starting values for the regression coefficients
#'   \eqn{\beta}.
#' @param data an optional data frame.
#' @param subset,na.action as in [stats::lm()].
#' @param knots the knots of \eqn{H_0(t)}: a number (the number of internal
#'   knots, placed at the quantiles of the observed event times) or a vector
#'   with the positions of all knots, including the boundaries (the first knot
#'   must not be larger than the smallest time, the last not smaller than the
#'   largest time). By default, `max(1, min(floor(n.events / 30), 3))` internal
#'   knots are used.
#' @param wtau an optional weighting function \eqn{w(\tau)}, e.g. [wtrunc()].
#' @param control a list created by [Qcoxph.control()].
#' @param ... additional arguments of `wtau`.
#'
#' @details
#' The model is \eqn{h(t \mid x) = h_0(t)\exp(x^\top\beta)}, with
#' \eqn{H_0(t) = \sum_l \exp(\gamma_l)\,b_l(t)}, where \eqn{b_l} is a piecewise
#' linear basis. All quantities of the estimating equations have closed forms
#' (Appendix D2 of the paper). The coefficients have the same interpretation
#' as in [survival::coxph()]; the estimates are typically very close, but more
#' robust. Starting points are obtained from the partial likelihood and the
#' Breslow estimator of \eqn{H_0}.
#'
#' Methods are available to compute the fitted survival curves ([survfit.Qcoxph()],
#' which returns an object of class `"survfit"`), predictions
#' ([predict.Qcoxph()]), residuals ([residuals.Qcoxph()]) and the baseline
#' cumulative hazard ([basehaz.Qcoxph()]). All of them are based on the fitted
#' parametric model.
#'
#' @return An object of class `c("Qcoxph", "Qest")`, a list with items
#' \item{coefficients, var}{the estimated \eqn{\beta} and their covariance matrix.}
#' \item{obj.function}{the loss at the solution.}
#' \item{iter, converged}{number of iterations and convergence.}
#' \item{linear.predictors}{\eqn{x_i^\top\beta - \bar x^\top\beta}.}
#' \item{residuals}{martingale residuals.}
#' \item{means, n, nevent}{means of the covariates, number of observations and of events.}
#' \item{concordance}{Harrell's concordance index (from [survival::concordancefit()]).}
#' \item{internal}{the knots, the coefficients \eqn{\gamma} of \eqn{H_0}, the full
#'   covariance matrix (`covar`), \eqn{F(y_i \mid x_i)} (`CDF`), \eqn{f(y_i \mid x_i)}
#'   (`PDF`), and other quantities used by the methods.}
#' and `terms`, `assign`, `formula`, `call`, `y`, `x`, `model`, `xlevels`,
#' `na.action`.
#' @references Sottile, G. and Frumento, P. (2022). Robust estimation and
#'   regression with parametric quantile functions. *Computational Statistics
#'   and Data Analysis*, 171, 107471. \doi{10.1016/j.csda.2022.107471}
#' @seealso [survfit.Qcoxph()], [predict.Qcoxph()], [residuals.Qcoxph()],
#'   [basehaz.Qcoxph()], [Qest()]
#' @examples
#' library(survival)
#' set.seed(1)
#' n <- 200
#' x <- runif(n, 0, 3)
#' t <- rweibull(n, shape = 2, scale = (1 / exp(2 + 2 * x))^(1 / 2))
#' cens <- runif(n, 0, 1)
#' dat <- data.frame(y = pmin(t, cens), d = as.numeric(t <= cens), x = x)
#'
#' m1 <- coxph(Surv(y, d) ~ x, data = dat)
#' m2 <- Qcoxph(Surv(y, d) ~ x, data = dat)
#' summary(m2)
#' cbind(coxph = coef(m1), Qcoxph = coef(m2))
#' sf <- survfit(m2, newdata = data.frame(x = c(0.5, 1.5)))
#' plot(sf, col = 1:2)
#' @export
Qcoxph <- function(formula, weights, start, data, subset, na.action, knots, wtau = NULL,
                   control = Qcoxph.control(), ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  na.act <- attr(mf, "na.action")
  if (!is.null(attr(Terms, "specials")$strata) || length(grep("strata\\(|cluster\\(|tt\\(", deparse(formula))))
    stop("strata(), cluster() and tt() terms are not supported", call. = FALSE)
  Y <- stats::model.response(mf)
  if (!inherits(Y, "Surv")) stop("the response must be a 'Surv' object", call. = FALSE)
  w <- check_weights(stats::model.weights(mf), nrow(mf))
  if (any(w == 0)) {
    warning("observations with zero weight are dropped", call. = FALSE)
    keep <- w > 0
    mf <- mf[keep, , drop = FALSE]; Y <- Y[keep, ]; w <- w[keep]
  }
  w <- w / mean(w)
  if ((n <- nrow(mf)) == 0) stop("zero non-NA cases", call. = FALSE)
  if (!is.null(stats::model.offset(mf))) stop("offsets are not supported", call. = FALSE)

  X <- stats::model.matrix(Terms, mf)
  Xatt <- attributes(X)
  xdrop <- Xatt$assign %in% 0
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  if (ncol(X) == 0) stop("the model must contain at least one covariate", call. = FALSE)
  nms <- colnames(X)

  resp <- get_response(Y)
  y <- resp$y; z <- resp$z; d <- resp$d; type <- resp$type
  if (any(y <= 0)) stop("times must be positive", call. = FALSE)

  dots <- list(...)
  if (!is.null(wtau)) {
    if (!is.function(wtau)) stop("'wtau' must be a function", call. = FALSE)
    split_dots(dots, wtau)
  } else if (length(dots)) stop("arguments in '...' are only allowed with 'wtau'", call. = FALSE)
  control <- do.call(Qest.control, control)
  tau <- cox_tau(wtau, dots)

  # knots ##################################################################
  ry <- range(y)
  ev <- y[d == 1]
  if (missing(knots) || is.null(knots)) {
    k <- max(1, min(floor(sum(d) / 30), 3))
    knots <- c(ry[1] - 1e-6, stats::quantile(ev, (1:k) / (k + 1), names = FALSE), ry[2] + 1e-6)
  } else if (length(knots) == 1) {
    if (knots < 0 || knots != round(knots)) stop("'knots' must be a non-negative integer or a vector", call. = FALSE)
    knots <- c(ry[1] - 1e-6, if (knots > 0) stats::quantile(ev, (1:knots) / (knots + 1), names = FALSE), ry[2] + 1e-6)
  } else {
    knots <- sort(unique(knots))
    if (any(!is.finite(knots))) stop("all knots must be finite", call. = FALSE)
    if (knots[1] > min(y) || knots[length(knots)] < max(y))
      stop("all times must be between the first and the last knot", call. = FALSE)
    test <- cut(ev, knots, include.lowest = TRUE)
    if (any(table(test) == 0)) stop("there must be some event between any two knots", call. = FALSE)
  }
  if (any(table(cut(ev, knots, include.lowest = TRUE)) == 0))
    stop("there must be some event between any two knots: use fewer knots", call. = FALSE)
  z <- pmax(z, knots[1])   # for z < knots[1], H0(z) would be negative

  # scaling, singularities, starting points ##############################
  X0 <- X; z0 <- z; y0 <- y; knots0 <- knots
  sv <- scalevars.Qcoxph(X, z, y, knots)
  if (any(!is.finite(sv$X))) stop("some covariates are constant", call. = FALSE)
  ok <- independent_cols(sv$X, tol = 1e-10)
  if (any(!ok)) warning("some variables were dropped because of collinearity: ", paste(nms[!ok], collapse = ", "),
                        call. = FALSE)
  sv$Stats$mX <- sv$Stats$mX[ok]; sv$Stats$sX <- sv$Stats$sX[ok]
  Xs <- sv$X[, ok, drop = FALSE]; zs <- sv$z; ys <- sv$y; ks <- sv$knots
  rank <- sum(ok)

  theta <- cox_start(Xs, zs, ys, d, w, ks)
  if (!missing(start) && !is.null(start)) {
    if (length(start) != ncol(X)) stop("length(start) does not match the number of predictors", call. = FALSE)
    if (anyNA(start)) stop("NAs are not allowed in 'start'", call. = FALSE)
    theta[1:rank] <- start[ok] * sv$Stats$sX
  }

  ee <- if (type == "ct") QCox.ee.ct else QCox.ee.c
  eefun <- function(th, J, zz = zs, yy = ys, XX = Xs, kk = ks) {
    G <- suppressWarnings(try(ee(th, 0, zz, yy, d, XX, w, kk, tau, J = J), silent = TRUE))
    if (inherits(G, "try-error") || any(!is.finite(G$g)) || (J && any(!is.finite(G$J)))) return(NULL)
    if (J) G$J <- t(G$J)   # the Jacobian is returned transposed (as in the families)
    G
  }
  fit <- qest_solve(theta, eefun, n, control)

  # back to the original scale ##########################################
  theta <- descalecoef.Qcoxph(fit$theta, sv$Stats)
  Xok <- X0[, ok, drop = FALSE]
  G <- eefun(theta, TRUE, zz = z0, yy = y0, XX = Xok, kk = knots0)
  if (is.null(G)) stop("the estimating equations cannot be evaluated at the solution", call. = FALSE)
  covar <- tryCatch(sandwich(G$J, G$g.i, w), error = function(e) NULL)
  if (is.null(covar)) stop("could not fit the model (singular Jacobian): ",
                           "try changing the number or the position of the knots", call. = FALSE)
  L <- coxLoss(theta, z0, y0, d, Xok, w, knots0, tau, type, G$Fy, G$Fz)
  attr(L, "npar") <- length(theta)
  if (type != "u") attr(L, "message") <- "Not the function being minimized"

  # output ##############################################################
  gnames <- paste0("log(slope", seq_len(length(theta) - rank), ")")
  pn <- c(nms[ok], gnames)
  dimnames(covar) <- list(pn, pn)
  beta <- setNames(rep(NA_real_, ncol(X)), nms); beta[ok] <- theta[1:rank]
  gamma <- theta[-(1:rank)]; names(gamma) <- gnames
  var <- matrix(0, ncol(X), ncol(X), dimnames = list(nms, nms))
  var[ok, ok] <- covar[1:rank, 1:rank, drop = FALSE]
  means <- colMeans(X0)
  eta <- drop(Xok %*% theta[1:rank])
  H0y <- c(plfcox(y0, knots0) %*% exp(gamma))
  H0z <- if (type == "ct") c(plfcox(z0, knots0) %*% exp(gamma)) else 0
  resid <- d - (H0y - H0z) * exp(eta)   # martingale residuals
  lp <- eta - sum(theta[1:rank] * means[ok])
  names(lp) <- names(resid) <- rownames(mf)
  conc <- cox_concordance(Y, lp)
  Jm <- G$J / n; dimnames(Jm) <- list(pn, pn)
  eeg <- G$g / n; names(eeg) <- pn

  internal <- list(z = z0, y = y0, d = d, w = w, type = type, tau = tau, knots = knots0, gamma = gamma,
                   var_gamma = covar[-(1:rank), -(1:rank), drop = FALSE], covar = covar, ee = eeg, jacobian = Jm,
                   CDF = G$Fy, PDF = G$fy, converged = fit$converged, ok = ok, wtau = wtau, wtauoptions = dots)
  out <- list(coefficients = beta, var = var, obj.function = L, iter = fit$n.it, converged = fit$converged,
              linear.predictors = lp, residuals = resid, means = means, n = n, nevent = sum(d),
              terms = Terms, assign = cox_assign(attr(X, "assign"), Terms), concordance = conc,
              model = mf, x = X0, y = Y, formula = formula, call = cl, xlevels = stats::.getXlevels(Terms, mf),
              na.action = na.act, internal = internal)
  class(out) <- c("Qcoxph", "Qest")
  if (!fit$converged) warning("the algorithm did not converge", call. = FALSE)
  out
}

# Functions of tau used by the closed forms of Qcoxph (w, W, its integral, and the
# corresponding r(tau*) functions of tau* = log(-log(1 - tau))).
cox_tau <- function(wtau, opt) {
  tauL <- 1e-10; tauR <- 1 - 1e-10
  taustarL <- log(-log(1 - tauL)); taustarR <- log(-log(1 - tauR))
  if (is.null(wtau)) {
    wfun <- function(tau, ...) tau * 0 + 1
    Wfun <- function(tau) tau
    iWfun <- function(tau) 0.5 * tau^2
    w1fun <- function(tau) 0 * tau
    rfun <- function(taustar, ...) exp(taustar)   # -log(1 - tau)
    Rfun <- function(taustar) {
      etaustar <- exp(taustar); tau1 <- exp(-etaustar)
      pmax0(1 - tau1 * (1 + etaustar))            # tau + (1 - tau) log(1 - tau)
    }
    iRfun <- function(taustar) {
      etaustar <- exp(taustar); tau1 <- exp(-etaustar); tau <- 1 - tau1
      pmax0(0.75 * tau^2 - 0.5 * (tau - (tau1^2) * etaustar))
    }
    r1fun <- function(tau, ...) 1 / (1 - tau)
  } else {
    # w(tau) and r(tau) are evaluated exactly; W, R and their integrals are
    # interpolated (monotone splines); w1 = w' and r1 = r' by finite differences,
    # taking the side with the smallest absolute value (for discontinuous w).
    ntau <- 4999
    tau <- seq(tauL, tauR, length.out = ntau); dtau <- (tau[2] - tau[1]) * 0.5
    taustar <- seq(taustarL, taustarR, length.out = ntau); dtaustar <- (taustar[2] - taustar[1]) * 0.5
    h <- exp(taustar - exp(taustar))
    wfun <- function(tau, ...) do.call(wtau, c(list(tau), opt))
    wt <- wfun(tau)
    if (anyNA(wt) || any(wt < 0) || any(!is.finite(c(wt, wfun(0), wfun(1)))))
      stop("'wtau(tau)' must return a non-negative, finite value for 0 <= tau <= 1", call. = FALSE)
    Wtau <- cumsum(wt[1:(ntau - 1)] + wt[2:ntau]) * dtau
    iWtau <- cumsum(Wtau[1:(ntau - 2)] + Wtau[2:(ntau - 1)]) * dtau
    w1tau <- (wt[-1] - wt[-ntau]) / (2 * dtau)
    Wfun <- stats::splinefun(c(0, tau[-1]), c(0, Wtau), method = "hyman")
    iWfun <- stats::splinefun(c(0, tau[-(1:2)]), c(0, iWtau), method = "hyman")
    w1fun <- stats::approxfun(tau[-1], minabs(w1tau, c(w1tau[-1], Inf)), rule = 2)
    rfun <- function(taustar, ...) {
      etaustar <- exp(taustar); tau1 <- exp(-etaustar)
      etaustar * wfun(1 - tau1)
    }
    rtauh <- rfun(taustar) * h
    Rtau <- cumsum(rtauh[1:(ntau - 1)] + rtauh[2:ntau]) * dtaustar
    Rtauh <- Rtau * h[-1]
    iRtau <- cumsum(Rtauh[1:(ntau - 2)] + Rtauh[2:(ntau - 1)]) * dtaustar
    Rfun <- stats::splinefun(c(taustarL, taustar[-1]), c(0, Rtau), method = "hyman")
    iRfun <- stats::splinefun(c(taustarL, taustar[-(1:2)]), c(0, iRtau), method = "hyman")
    r1fun <- function(tau, ...) { tau <- pmin(tau, 1 - 1e-10); wfun(tau) / (1 - tau) - log(1 - tau) * w1fun(tau) }
  }
  list(wfun = wfun, Wfun = Wfun, iWfun = iWfun, w1fun = w1fun,
       rfun = rfun, Rfun = Rfun, iRfun = iRfun, r1fun = r1fun,
       tauL = tauL, tauR = tauR, taustarL = taustarL, taustarR = taustarR, opt = list())
}

minabs <- function(x1, x2) {
  i <- which(abs(x1) < abs(x2))
  out <- x2; out[i] <- x1[i]
  out
}

# Starting points of Qcoxph: beta from the partial likelihood (Breslow ties, left
# truncation allowed), H0 from the Breslow estimator, fitted by a piecewise linear
# function with positive slopes.
cox_start <- function(X, z, y, d, w, knots) {
  p <- ncol(X)
  ev <- which(d == 1)
  te <- sort(unique(y[ev]))
  oy <- order(y); ys <- y[oy]
  zf <- is.finite(z); oz <- order(z); zsrt <- z[oz]
  # sums over the risk set {i: z_i < t <= y_i}, for each event time t
  risksum <- function(v) {   # v: n-vector or n-column matrix
    v <- as.matrix(v)
    cy <- apply(v[oy, , drop = FALSE], 2, function(u) rev(cumsum(rev(u))))
    cy <- rbind(as.matrix(cy), 0)
    iy <- findInterval(te, ys, left.open = TRUE) + 1   # first index with y >= t
    out <- cy[iy, , drop = FALSE]
    if (any(zf)) {
      cz <- apply(v[oz, , drop = FALSE], 2, function(u) rev(cumsum(rev(u))))
      cz <- rbind(as.matrix(cz), 0)
      iz <- findInterval(te, zsrt, left.open = TRUE) + 1    # first index with z >= t
      out <- out - cz[iz, , drop = FALSE]
    }
    out
  }
  de <- as.numeric(tapply(w[ev], factor(y[ev], levels = te), sum))
  sxe <- rowsum(w[ev] * X[ev, , drop = FALSE], factor(y[ev], levels = te), reorder = TRUE)
  loglik <- function(b) {
    r <- w * exp(c(X %*% b))
    s0 <- risksum(r)[, 1]
    sum(sxe %*% b) - sum(de * log(s0))
  }
  beta <- rep(0, p)
  ll <- loglik(beta)
  for (it in 1:30) {
    r <- w * exp(c(X %*% beta))
    s0 <- risksum(r)[, 1]
    s1 <- risksum(r * X)
    XX <- tensorX(X)
    s2 <- risksum(r * XX)
    xbar <- s1 / s0
    grad <- colSums(sxe) - colSums(de * xbar)
    H <- matrix(0, p, p)
    h2 <- colSums(de * s2 / s0)
    H <- sym_from_lower(h2, p) - crossprod(xbar * sqrt(de))
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step)) break
    a <- 1
    repeat {
      nb <- beta + a * step
      nll <- loglik(nb)
      if (is.finite(nll) && nll >= ll - 1e-10) break
      a <- a / 2
      if (a < 1e-6) { nb <- beta; nll <- ll; break }
    }
    conv <- max(abs(nb - beta)) < 1e-8
    beta <- nb; ll <- nll
    if (conv) break
  }
  r <- w * exp(c(X %*% beta))
  s0 <- risksum(r)[, 1]
  H0 <- cumsum(de / s0)
  b0 <- plfcox(te, knots)
  slope <- stats::lm.fit(b0, H0)$coefficients
  slope[!is.finite(slope)] <- 0
  slope <- adjust.coef(slope)
  if (any(slope <= 0)) slope[slope <= 0] <- max(c(H0[length(H0)] / diff(range(knots)), 1e-6))
  c(beta, log(slope))
}

cox_concordance <- function(Y, lp) {
  cc <- tryCatch(survival::concordancefit(Y, lp, reverse = TRUE, timefix = FALSE), error = function(e) NULL)
  if (is.null(cc)) return(NULL)
  c(cc$count, concordance = cc$concordance, std = sqrt(cc$var))
}

cox_assign <- function(asgn, Terms) {
  lab <- attr(Terms, "term.labels")
  out <- lapply(seq_along(lab), function(j) which(asgn == j))
  names(out) <- lab
  out[lengths(out) > 0]
}
