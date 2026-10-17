#' Q-Estimation of Linear Regression Models
#'
#' Fits a linear model with normal homoskedastic errors,
#' \eqn{Q(\tau \mid x) = x^\top\beta + \sigma\Phi^{-1}(\tau)}, with the
#' Q-estimator of Sottile and Frumento (2022). The estimating equations and the
#' Hessian have closed forms (Appendix D1 of the paper), and the estimator is
#' fast and robust to outliers in the response.
#'
#' @param formula,data,subset,weights,na.action,contrasts as in [stats::lm()].
#' @param start optional starting values for the regression coefficients
#'   (excluding \eqn{\sigma}).
#' @param wtau an optional weighting function \eqn{w(\tau)}, e.g. [wtrunc()].
#' @param control a list created by [Qest.control()].
#' @param ... additional arguments of `wtau`.
#' @param y,X,w for `Qlm.fit`: the response, the design matrix and the weights.
#'
#' @details
#' `Qlm(formula)` estimates the same model as `Qest(formula, Q = Qnorm)`, with
#' a much richer output. Censored and truncated data are not allowed (use
#' [Qest()] with `Q = Qnorm`). Standard errors are computed with the sandwich
#' formula, and the tests in [summary.Qlm()] are based on the asymptotic normal
#' distribution.
#'
#' @return `Qlm` returns an object of class `c("Qlm", "Qest")` with items
#' \item{coefficients, std.errs, covar}{the estimated regression coefficients,
#'   their standard errors and covariance matrix.}
#' \item{dispersion}{the estimated \eqn{\sigma^2}.}
#' \item{sigma}{the estimated \eqn{\sigma}, with attribute `"std.err"` (the
#'   standard error of \eqn{\log\sigma}).}
#' \item{residuals, fitted.values}{residuals and fitted values.}
#' \item{CDF}{the fitted \eqn{F(y_i \mid x_i)}.}
#' \item{obj.function}{the minimized loss.}
#' \item{gradient, hessian}{the estimating equations and the Hessian of the loss
#'   at the solution (including \eqn{\log\sigma}).}
#' \item{covar.full}{the covariance matrix of \eqn{(\log\sigma, \beta)}.}
#' and `rank`, `df.residual`, `weights`, `converged`, `n.it`, `call`, `terms`,
#' `model`, `xlevels`, `contrasts`, `na.action`. `Qlm.fit` returns a similar list.
#' @references Sottile, G. and Frumento, P. (2022). Robust estimation and
#'   regression with parametric quantile functions. *Computational Statistics
#'   and Data Analysis*, 171, 107471. \doi{10.1016/j.csda.2022.107471}
#' @seealso [summary.Qlm()], [predict.Qlm()], [Qest()]
#' @examples
#' set.seed(1234)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- runif(n, 0, 3)
#' y <- rnorm(n, 4 + x1 + 2 * x2, 1)
#' y[1:5] <- y[1:5] + 20      # outliers
#'
#' m1 <- Qlm(y ~ x1 + x2)
#' summary(m1)
#' m2 <- Qlm(y ~ x1 + x2, wtau = wtrunc, delta.left = 0.05, delta.right = 0.05)
#' cbind(lm = coef(lm(y ~ x1 + x2)), Qlm = coef(m1), Qlm.trunc = coef(m2))
#' @export
Qlm <- function(formula, data, subset, weights, na.action, start = NULL, contrasts = NULL,
                wtau = NULL, control = Qest.control(), ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (stats::is.empty.model(mt)) stop("the model has no covariates and no intercept", call. = FALSE)
  if (!is.null(stats::model.offset(mf))) stop("offsets are not supported", call. = FALSE)
  y <- stats::model.response(mf, "numeric")
  if (inherits(y, "Surv") || is.matrix(y)) stop("the response must be a numeric vector (use Qest with Q = Qnorm for censored data)", call. = FALSE)
  w <- check_weights(stats::model.weights(mf), nrow(mf))
  na.act <- attr(mf, "na.action")
  if (any(w == 0)) {
    warning("observations with zero weight are dropped", call. = FALSE)
    keep <- w > 0
    mf <- mf[keep, , drop = FALSE]; y <- y[keep]; w <- w[keep]
  }
  if (nrow(mf) == 0) stop("zero non-NA cases", call. = FALSE)
  x <- stats::model.matrix(mt, mf, contrasts)
  z <- Qlm.fit(y, x, w, start, wtau, control, ...)
  z$na.action <- na.act
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- stats::.getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  class(z) <- c("Qlm", "Qest")
  if (!z$converged) warning("the algorithm did not converge", call. = FALSE)
  z
}

#' @rdname Qlm
#' @export
Qlm.fit <- function(y, X, w = rep(1, nobs), start = NULL, wtau = NULL, control = Qest.control(), ...) {
  control <- do.call(Qest.control, control)
  X <- as.matrix(X)
  nobs <- nrow(X); nvars <- ncol(X)
  if (length(y) != nobs) stop("'y' and 'X' have different numbers of observations", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(nvars))
  w <- check_weights(w, nobs); w <- w / mean(w)
  if (!is.null(wtau)) {
    if (!is.function(wtau)) stop("'wtau' must be a function", call. = FALSE)
    split_dots(list(...), wtau)
  } else if (length(list(...))) stop("arguments in '...' are only allowed with 'wtau'", call. = FALSE)
  if (is.null(wtau)) wtau <- function(tau, ...) rep(1, length(tau))

  ok <- independent_cols(X)
  rank <- sum(ok)
  xnames <- colnames(X)
  if (any(!ok)) warning("some variables were dropped because of collinearity: ", paste(xnames[!ok], collapse = ", "),
                        call. = FALSE)
  X <- X[, ok, drop = FALSE]
  X0 <- X; y0 <- y
  sv <- scalevars.Qlm(X, y)
  X <- sv$X; y <- sv$y; ss <- sv$Stats

  # starting points: median regression and a robust estimate of sigma
  m0 <- rq_fnb(X, y, tau = 0.5)
  s0 <- diff(stats::quantile(m0$residuals, c(0.3, 0.7), names = FALSE)) / diff(stats::qnorm(c(0.3, 0.7)))
  theta0 <- c(log(max(s0, 1e-3)), m0$coefficients)
  if (!is.null(start)) {
    if (length(start) != nvars) stop("length(start) does not match the number of predictors", call. = FALSE)
    if (anyNA(start)) stop("NAs are not allowed in 'start'", call. = FALSE)
    st <- start[ok]
    st[ss$int] <- st[ss$int] - ss$my * ss$is.int
    st <- st * 10 / (ss$My - ss$my)
    st[ss$vars] <- (st * ss$sX)[ss$vars]
    st[ss$int] <- st[ss$int] + sum((st * ss$mX / ss$sX)[ss$vars])
    theta0[-1] <- st
  }

  bfun <- Qlm.bfun(wtau, ...)
  eefun <- function(th, J, yy = y, XX = X) {
    L <- qlmLoss(th, yy, XX, w, bfun)
    G <- Qlm.ee.u(th, XX, w, bfun, L, J = J)
    if (any(!is.finite(G$g)) || (J && any(!is.finite(G$H)))) return(NULL)
    list(g = G$g, g.i = G$g.i, J = G$H, loss = L)
  }
  fit <- qest_solve(theta0, eefun, nobs, control)

  theta1 <- descalecoef.Qlm(fit$theta, ss)
  G <- eefun(theta1, TRUE, yy = y0, XX = X0)
  H <- G$J
  covar <- tryCatch(sandwich(H, G$g.i, w), error = function(e) NULL)
  if (is.null(covar)) stop("the Hessian is singular", call. = FALSE)
  pn <- c("log(sigma)", xnames[ok])
  dimnames(covar) <- dimnames(H) <- list(pn, pn)

  coefficients <- std.errs <- setNames(rep(NA_real_, nvars), xnames)
  coefficients[ok] <- theta1[-1]
  std.errs[ok] <- sqrt(diag(covar)[-1])
  fitted.values <- drop(X0 %*% theta1[-1])
  residuals <- drop(y0 - fitted.values)
  names(fitted.values) <- names(residuals) <- names(y0)
  sigma <- exp(theta1[1]); names(sigma) <- NULL
  attr(sigma, "std.err") <- sqrt(covar[1, 1])
  gradient <- G$g / nobs; names(gradient) <- pn
  CDF <- G$loss$p.i; names(CDF) <- names(y0)

  list(coefficients = coefficients, std.errs = std.errs, covar = covar[-1, -1, drop = FALSE],
       dispersion = sigma^2, sigma = sigma, residuals = residuals, fitted.values = fitted.values,
       rank = rank, df.residual = nobs - rank, weights = w, obj.function = G$loss$Loss,
       gradient = gradient, hessian = H, covar.full = covar, CDF = CDF,
       converged = fit$converged, n.it = fit$n.it, control = control, wtau = wtau, wtauoptions = list(...))
}
