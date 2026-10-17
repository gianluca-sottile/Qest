#' Family Objects for Qest
#'
#' Built-in parametric models for [Qest()], with analytical (or nearly
#' analytical) estimating equations: much faster and more accurate than the
#' equivalent user-defined quantile functions.
#'
#' @param offset an optional vector of offsets of the log-rate (`Qpois` only).
#'   An `offset()` term in the formula can be used instead.
#' @param min logical (`Qunif` only): if `TRUE`, fit a \eqn{U(a, b)} model; if
#'   `FALSE`, a \eqn{U(0, b)} model.
#'
#' @details
#' A family can be passed to [Qest()] by name (`"Qpois"`), as a function
#' (`Qpois`) or as a call (`Qpois()`); the latter allows to pass arguments.
#'
#' * `Qnorm`: normal homoskedastic model, \eqn{Q(\tau \mid x) = x^\top\beta + \sigma\Phi^{-1}(\tau)};
#'   parameters `log(sigma)`, `beta`. With complete data, [Qlm()] fits the
#'   same model with a richer output.
#' * `Qgamma`: Gamma model with shape \eqn{a} and log-scale \eqn{x^\top\beta};
#'   parameters `log(shape)`, `beta`.
#' * `Qpois`: Poisson model with log-rate \eqn{x^\top\beta}. To obtain a continuous
#'   quantile function, `qpois` is replaced by the inverse with respect to \eqn{y}
#'   of the regularized upper incomplete gamma function \eqn{\Gamma(y, \lambda)/\Gamma(y)}
#'   (Ilienko, 2013). The model should be applied to a jittered response,
#'   `y + runif(n)`.
#' * `Qunif`: Uniform model \eqn{U(a(x), b(x))} with \eqn{a(x) = x^\top\alpha} and
#'   \eqn{b(x) = x^\top\beta}, or \eqn{U(0, b(x))} if `min = FALSE`.
#'
#' `Qnorm` and `Qgamma` can be used with censored and truncated data; `Qpois`
#' and `Qunif` cannot. `Qgamma` and `Qpois` require an intercept.
#'
#' @return An object of class `"Qfamily"`.
#' @references Ilienko, A. (2013). Continuous counterparts of Poisson and binomial
#'   distributions and their properties. *Annales Univ. Sci. Budapest., Sect. Comp.*,
#'   39, 137-147.
#'
#'   Sottile, G. and Frumento, P. (2022). Robust estimation and regression with
#'   parametric quantile functions. *Computational Statistics and Data Analysis*,
#'   171, 107471. \doi{10.1016/j.csda.2022.107471}
#' @seealso [Qest()]
#' @examples
#' set.seed(1)
#' n <- 250
#' x <- runif(n)
#' eta <- 1 + 2 * x
#'
#' m1 <- Qest(y ~ x, Qnorm, data = data.frame(y = rnorm(n, eta, exp(1)), x = x))
#' m2 <- Qest(y ~ x, Qgamma, data = data.frame(y = rgamma(n, shape = exp(1), scale = exp(eta)), x = x))
#' m3 <- Qest(y ~ x, Qpois, data = data.frame(y = rpois(n, exp(eta)) + runif(n), x = x))
#' m4 <- Qest(y ~ x, Qunif(min = TRUE), data = data.frame(y = runif(n, 0, eta), x = x))
#' m5 <- Qest(y ~ x, Qunif(min = FALSE), data = data.frame(y = runif(n, 0, eta), x = x))
#' coef(m2)
#' @name Qfamily
NULL

#' @export
print.Qfamily <- function(x, ...) {
  nm <- x$family$family
  nm <- switch(nm, gaussian = "Qnorm", Gamma = "Qgamma", poisson = "Qpois", uniform = "Qunif", nm)
  cat("Qfamily:", nm, if (nm == "Qunif") paste0("(min = ", x$min, ")"), "\n")
  invisible(x)
}

# Fit with a Qfamily.
qest_family_fit <- function(Qf, z, y, d, x, w, type, start, ntau, wtau, wtauoptions, control, offset = NULL) {
  fam <- Qf$family$family
  if (fam %in% c("uniform", "poisson") && type != "u")
    stop("censored and truncated data are not supported by ", if (fam == "poisson") "Qpois" else "Qunif",
         call. = FALSE)
  nobs <- length(y)
  if (fam == "poisson") {
    if (!is.null(Qf$offset)) {
      if (length(Qf$offset) != nobs)
        stop(gettextf("the offset has length %d, but there are %d observations", length(Qf$offset), nobs),
             call. = FALSE)
      offset <- if (is.null(offset)) Qf$offset else offset + Qf$offset
    }
    if (is.null(offset)) offset <- rep.int(0, nobs)
  }
  tau <- buildTau(ntau, wtau, nobs, wtauoptions)
  if (fam == "poisson") tau <- tau.pois(tau)

  nms <- colnames(x)
  ok <- independent_cols(x)
  x <- x[, ok, drop = FALSE]
  if (any(!ok)) warning("some variables were dropped because of collinearity: ",
                        paste(nms[!ok], collapse = ", "), call. = FALSE)
  attr(x, "offset") <- offset
  attr(x, "min") <- Qf$min
  temp <- Qf$scale(x, y, z)
  colnames(temp$Xc) <- colnames(temp$Stats$X) <- nms[ok]
  Stats <- temp$Stats
  x <- temp$Xc; y <- temp$yc; z <- temp$zc
  attr(x, "offset") <- offset
  attr(x, "min") <- Qf$min

  theta <- if (is.null(start)) Qf$initialize(x, z, y, d, w, Qf, tau = tau, ok = ok, Stats = Stats)
           else Qf$initialize(x, z, y, d, w, Qf, start = start, tau = tau, ok = ok, Stats = Stats)
  pnames <- names(theta)

  Q0 <- Qf$Q
  attr(Q0, "glmFamily") <- Qf$family
  if (!is.null(Qf$bfun)) tau$bfun <- Qf$bfun(wtau, wtauoptions)
  attr(Q0, "X") <- x
  class(Q0) <- class(Qf)

  ee <- switch(type,
               u = switch(fam, Gamma = QestGamma.ee.u, poisson = QestPois.ee.u,
                          uniform = QestUnif.ee.u, gaussian = QestNorm.ee.u),
               c = switch(fam, Gamma = QestGamma.ee.c, gaussian = QestNorm.ee.c),
               ct = switch(fam, Gamma = QestGamma.ee.ct, gaussian = QestNorm.ee.ct))
  eefun <- function(th, J, Q = Q0, yy = y, zz = z) {
    G <- suppressWarnings(try(ee(th, 1, zz, yy, d, Q, w, NULL, tau, J = J), silent = TRUE))
    if (inherits(G, "try-error") || any(!is.finite(G$g)) || (J && any(!is.finite(G$J)))) return(NULL)
    # with censored or truncated data, the Jacobian of the families (not symmetric) was
    # returned transposed in Qest <= 1.0.2
    if (J && type != "u") G$J <- t(G$J)
    G
  }
  fit <- qest_solve(theta, eefun, nobs, control)

  # back to the original scale
  theta <- Qf$descale(fit$theta, Stats)
  Q1 <- Q0; attr(Q1, "X") <- Stats$X
  G <- eefun(theta, TRUE, Q = Q1, yy = Stats$y, zz = Stats$z)
  if (is.null(G)) stop("the estimating equations cannot be evaluated at the solution", call. = FALSE)
  if (fam == "poisson") { Fy <- G$Fy; Fz <- NULL }
  else {
    Fy <- Qf$fix.Fy(list(Fy = G$Fy), theta, tau, Stats$y, Stats$X, Q1)
    Fz <- if (type == "ct") Qf$fix.Fy(list(Fy = G$Fz), theta, tau, Stats$z, Stats$X, Q1)
  }
  loss <- family_loss(w, d, tau, type, Fy, Fz)
  attr(loss, "npar") <- length(theta)
  if (type != "u") attr(loss, "message") <- "Not the function being minimized"

  covar <- tryCatch(sandwich(G$J, G$g.i, w), error = function(e) NULL)
  if (is.null(covar)) {
    warning("the Jacobian is singular: is the model identifiable?", call. = FALSE)
    covar <- matrix(NA_real_, length(theta), length(theta))
  }
  dimnames(covar) <- list(pnames, pnames)
  J <- G$J / nobs; dimnames(J) <- list(pnames, pnames)
  ee <- G$g / nobs; names(ee) <- pnames
  names(theta) <- pnames
  se <- sqrt(pmax(diag(covar), 0))

  # coefficients of the dropped variables are NA
  okfull <- switch(fam,
                   gaussian = , Gamma = c(TRUE, ok),
                   uniform = rep(ok, 1 + Qf$min),
                   ok)
  fullnames <- switch(fam,
                      gaussian = c("log(sigma)", nms), Gamma = c("log(shape)", nms),
                      uniform = if (Qf$min) c(paste0("min: ", nms), paste0("max: ", nms)) else paste0("max: ", nms),
                      nms)
  th <- se2 <- setNames(rep(NA_real_, length(okfull)), fullnames)
  th[okfull] <- theta; se2[okfull] <- se

  list(theta = th, se = se2, covar = covar, loss = loss, ee = ee, jacobian = J,
       CDF = Fy$Fy, PDF = G$fy, converged = fit$converged, n.it = fit$n.it, Q = Q1, tau = tau,
       family = Qf)
}

# Loss of the families, computed on the grid tau.
family_loss <- function(w, d, tau, type, Fy, Fz) {
  TAU <- tau$TAU
  Sy <- if (type != "u") 1 - Fy$Fy else 1
  deltay <- Fy$delta
  omegay <- (deltay <= 0)
  om <- if (type == "u") omegay else omegay + (1 - d) * omegay * (TAU - 1) / Sy
  if (type == "ct") {
    omegaz <- (Fz$delta <= 0)
    om <- om - omegaz - omegaz * (TAU - 1) / (1 - Fz$Fy) + TAU
  }
  out <- (deltay * (TAU - om)) %*% (tau$dtau * tau$wtau)
  2 * sum(w * out)
}
