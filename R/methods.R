# Methods for Qest and Qlm objects.

ci_table <- function(cf, se, parm, level) {
  pnames <- names(cf)
  if (missing(parm)) parm <- pnames
  else if (is.numeric(parm)) parm <- pnames[parm]
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ci[] <- cf[parm] + se[parm] %o% stats::qnorm(a)
  ci
}

#' @export
print.Qest <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(stats::coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' Summarizing Q-Estimators
#'
#' Summary methods for objects of class `"Qest"` (fitted by [Qest()]) and
#' `"Qlm"` (fitted by [Qlm()]).
#'
#' @param object an object of class `"Qest"` or `"Qlm"`.
#' @param covar logical: include the covariance matrix?
#' @param correlation logical: include the correlation matrix of the
#'   coefficients (`summary.Qlm`)?
#' @param x a summary object.
#' @param digits,signif.stars printing options.
#' @param ... not used.
#'
#' @details
#' The tests are based on the asymptotic normal distribution of the estimators,
#' with standard errors from the sandwich formula (Sottile and Frumento, 2022).
#'
#' @return `summary.Qest` returns an object of class `"summary.Qest"` with the
#'   table of coefficients (`coefficients`), the loss (`obj.function`), the
#'   numbers of observations and parameters, the convergence status and,
#'   optionally, the covariance matrix. `summary.Qlm` returns an object of class
#'   `"summary.Qlm"` with the table of coefficients, the estimate of
#'   \eqn{\sigma} and its standard error, the residuals and, optionally, the
#'   correlation matrix.
#' @seealso [Qest()], [Qlm()]
#' @examples
#' set.seed(1)
#' y <- rgamma(200, shape = 3, scale = 2)
#' fit <- Qest(y ~ 1, Q = Qgamma)
#' summary(fit, covar = TRUE)
#' @export
summary.Qest <- function(object, covar = FALSE, ...) {
  theta <- object$coefficients
  se <- object$std.errs
  coe <- cbind(theta, se, theta / se, 2 * stats::pnorm(-abs(theta / se)))
  dimnames(coe) <- list(names(theta), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  out <- list(coefficients = coe, obj.function = object$obj.function, n = length(object$internal$y),
              npar = sum(!is.na(theta)), iter = object$n.it, converged = object$converged,
              covar = if (covar) object$covar, call = object$call, type = object$internal$type)
  class(out) <- "summary.Qest"
  out
}

#' @rdname summary.Qest
#' @export
print.summary.Qest <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of observations:", x$n, "\n")
  cat("n. of parameters:", x$npar, "\n")
  cat("n. of iterations:", x$iter, if (!x$converged) "(not converged)", "\n\n")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  cat("\n", if (x$type == "u") "Minimized loss function: " else "Loss (not the function being minimized): ",
      format(x$obj.function, digits = digits + 3), "\n", sep = "")
  if (!is.null(x$covar)) {
    cat("\nCovariance matrix:\n")
    print.default(format(x$covar, digits = digits), print.gap = 2L, quote = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' @export
vcov.Qest <- function(object, ...) object$covar

#' @export
confint.Qest <- function(object, parm, level = 0.95, ...) {
  ci_table(stats::coef(object), object$std.errs, parm, level)
}

#' @export
nobs.Qest <- function(object, ...) length(object$internal$y)

#' Predictions from Q-Estimators
#'
#' Computes the fitted quantiles \eqn{Q(p \mid x)}, or the fitted distribution
#' function \eqn{F(y \mid x)} at the observed (or new) values of the response.
#'
#' @param object an object of class `"Qest"`.
#' @param newdata an optional data frame. For `type = "CDF"`, it must contain the
#'   response.
#' @param type `"quantile"` (default) or `"CDF"`.
#' @param p the quantile orders for `type = "quantile"`.
#' @param ... additional arguments of a user-defined `Q`, if they differ from
#'   those used to fit the model.
#' @return For `type = "quantile"`, a matrix with one row for each observation
#'   and one column for each value of `p`; for `type = "CDF"`, a vector.
#' @details For a user-defined `Q`, `newdata` is passed as the argument `data`
#'   of `Q`. For `type = "CDF"`, \eqn{F(y \mid x)} is computed by inverting
#'   \eqn{Q} (see [invQ()]). For censored data, the CDF is evaluated at the
#'   observed times.
#' @seealso [Qest()], [invQ()]
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = runif(200))
#' dat$y <- rgamma(200, shape = 2, scale = exp(1 + dat$x))
#' fit <- Qest(y ~ x, Q = Qgamma, data = dat)
#' predict(fit, newdata = data.frame(x = c(0, 1)), p = c(0.1, 0.5, 0.9))
#' head(predict(fit, type = "CDF"))
#' @export
predict.Qest <- function(object, newdata, type = c("quantile", "CDF"), p = c(0.1, 0.25, 0.5, 0.75, 0.9), ...) {
  type <- match.arg(type)
  if (type == "quantile" && any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
  if (missing(newdata)) {
    if (type == "CDF") return(object$CDF)
    newdata <- NULL
  }
  fam <- object$internal$family
  if (is.null(fam)) {
    Q <- object$internal$Q
    dots <- list(...)
    if (length(dots)) attr(Q, "Qoptions") <- utils::modifyList(attr(Q, "Qoptions"), dots)
    dat <- if (is.null(newdata)) object$data else newdata
    if (!is.data.frame(dat)) dat <- as.data.frame(dat)
    theta <- object$coefficients
    if (type == "quantile") {
      out <- callQ(Q, theta, matrix(p, nrow(dat), length(p), byrow = TRUE), dat)
      out <- matrix(out, nrow(dat), length(p), dimnames = list(rownames(dat), paste0("p", p)))
      return(out)
    }
    mf <- stats::model.frame(object$terms, dat, na.action = stats::na.pass)
    yy <- get_response(stats::model.response(mf))$y
    return(invQ(Q, theta, yy, dat))
  }
  # families
  famname <- fam$family$family
  tt <- stats::delete.response(object$terms)
  if (is.null(newdata)) X <- object$internal$x
  else {
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = object$xlevels)
    X <- stats::model.matrix(tt, mf)
  }
  th <- object$coefficients; th[is.na(th)] <- 0
  q <- ncol(X)
  pars <- switch(famname,
                 gaussian = list(sigma = exp(th[1]), mu = c(X %*% th[-1])),
                 Gamma = list(shape = exp(th[1]), scale = exp(c(X %*% th[-1]))),
                 poisson = list(lambda = exp(c(X %*% th) + (if (!is.null(fam$offset)) fam$offset else 0))),
                 uniform = list(a = if (fam$min) c(X %*% th[1:q]) else rep(0, nrow(X)),
                                b = c(X %*% th[(q * fam$min + 1):(q * fam$min + q)])))
  if (type == "quantile") {
    out <- sapply(p, function(pp) switch(famname,
      gaussian = stats::qnorm(pp, pars$mu, pars$sigma),
      Gamma = stats::qgamma(pp, shape = pars$shape, scale = pars$scale),
      poisson = qpois_cont(pp, pars$lambda),
      uniform = pars$a + pp * (pars$b - pars$a)))
    out <- matrix(out, nrow(X), length(p), dimnames = list(rownames(X), paste0("p", p)))
    return(out)
  }
  mf <- stats::model.frame(object$terms, newdata, na.action = stats::na.pass, xlev = object$xlevels)
  yy <- get_response(stats::model.response(mf))$y
  switch(famname,
         gaussian = stats::pnorm(yy, pars$mu, pars$sigma),
         Gamma = stats::pgamma(yy, shape = pars$shape, scale = pars$scale),
         poisson = ppoisC(yy, pars$lambda),
         uniform = pmin(pmax((yy - pars$a) / (pars$b - pars$a), 0), 1))
}

# Quantiles of the continuous Poisson distribution, by bisection.
qpois_cont <- function(p, lambda) {
  lo <- rep(0, length(lambda)); hi <- stats::qpois(p, lambda) + 2
  for (i in 1:60) {
    mid <- (lo + hi) / 2
    up <- ppoisC(mid, lambda) < p
    lo[up] <- mid[up]; hi[!up] <- mid[!up]
  }
  (lo + hi) / 2
}

#' Diagnostic Plot of Q-Estimators
#'
#' If the model is correctly specified, \eqn{F(y_i \mid x_i)} are a sample from
#' a \eqn{U(0, 1)} distribution. The plot compares their empirical distribution
#' (the Kaplan-Meier estimator, with censored or truncated data) with the
#' uniform distribution function.
#'
#' @param x an object of class `"Qest"` (or `"Qlm"`).
#' @param conf.int logical: add pointwise 95% confidence bands (Greenwood)?
#' @param ... additional graphical parameters passed to [graphics::plot()].
#' @return Invisibly, a data frame with the Kaplan-Meier estimate of the
#'   distribution of \eqn{F(y_i \mid x_i)}.
#' @seealso [Qest()]
#' @examples
#' set.seed(1)
#' y <- rgamma(300, shape = 2, scale = 3)
#' plot(Qest(y ~ 1, Q = Qgamma))   # good fit
#' plot(Qest(y ~ 1, Q = Qnorm))    # poor fit
#' @export
plot.Qest <- function(x, conf.int = TRUE, ...) {
  u <- x$CDF
  int <- x$internal
  if (is.null(int)) int <- list(d = rep(1, length(u)), w = rep(1, length(u)), type = "u")   # Qlm
  km <- km_unif(u, int$d, int$w, if (int$type == "ct") truncation_cdf(x) else NULL)
  args <- list(...)
  defaults <- list(x = c(0, 1), y = c(0, 1), type = "n", xlab = "U(0,1) quantiles",
                   ylab = expression(paste("empirical distribution of ", hat(F)(y[i] ~ "|" ~ x[i]))),
                   main = "Uniform probability plot")
  do.call(graphics::plot, utils::modifyList(defaults, args))
  graphics::abline(0, 1, lty = 2, col = "grey40")
  graphics::lines(km$u, km$F, type = "s", lwd = 1.5)
  if (conf.int) {
    graphics::lines(km$u, km$lower, type = "s", lty = 3)
    graphics::lines(km$u, km$upper, type = "s", lty = 3)
  }
  invisible(km)
}

# F(z_i | x_i) for truncated data (only used by plot.Qest).
truncation_cdf <- function(x) {
  int <- x$internal
  if (is.null(int$family)) return(invQ(int$Q, x$coefficients, int$z, x$data))
  famname <- int$family$family$family
  th <- x$coefficients; th[is.na(th)] <- 0
  X <- int$x
  if (famname == "gaussian") stats::pnorm(int$z, c(X %*% th[-1]), exp(th[1]))
  else stats::pgamma(int$z, shape = exp(th[1]), scale = exp(c(X %*% th[-1])))
}

# Kaplan-Meier estimator of the distribution of u (with censoring d and
# truncation uz), with Greenwood bands.
km_unif <- function(u, d, w, uz = NULL) {
  tt <- sort(unique(u[d == 1]))
  if (is.null(uz)) uz <- rep(-Inf, length(u))
  nr <- vapply(tt, function(t) sum(w[u >= t & uz < t]), 1)
  ne <- vapply(tt, function(t) sum(w[u == t & d == 1]), 1)
  S <- cumprod(1 - ne / nr)
  v <- cumsum(ne / (nr * pmax(nr - ne, 1e-10)))
  se <- S * sqrt(v)
  data.frame(u = c(0, tt), F = c(0, 1 - S), lower = c(0, pmax(1 - S - 1.96 * se, 0)),
             upper = c(0, pmin(1 - S + 1.96 * se, 1)))
}

#' Inverse of a Quantile Function
#'
#' Given a quantile function \eqn{Q(\tau \mid x, \theta)} and a value of
#' \eqn{\theta}, computes \eqn{F(y_i \mid x_i, \theta)} by bisection.
#'
#' @param Q a quantile function `Q(theta, tau, data)` (see [Qest()]).
#' @param theta the parameters.
#' @param y the response: one value for each row of `data`.
#' @param data a data frame.
#' @param n.it the number of iterations of the bisection algorithm (the accuracy
#'   is \eqn{2^{-n.it}}).
#' @return A vector with \eqn{F(y_i \mid x_i, \theta)}.
#' @seealso [Qest()]
#' @examples
#' Qlinmod <- function(theta, tau, data) {
#'   qnorm(tau, theta[2] + theta[3] * data$x, exp(theta[1]))
#' }
#' set.seed(1)
#' dat <- data.frame(x = runif(10))
#' U <- runif(10)
#' y <- Qlinmod(c(0, 1, 2), U, dat)
#' cbind(U, invQ(Qlinmod, c(0, 1, 2), y, dat))
#' @export
invQ <- function(Q, theta, y, data, n.it = 17) {
  n <- length(y)
  tau <- rep.int(0.5, n)
  for (i in 2:n.it) {
    omega <- sign(y - callQ(Q, theta, tau, data))
    tau <- tau + 1 / (2^i) * omega
  }
  tau
}

#### Qlm ##################################################################

#' @export
print.Qlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(stats::coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  cat("\nsigma:", format(x$sigma, digits = digits), "\n\n")
  invisible(x)
}

#' @rdname summary.Qest
#' @export
summary.Qlm <- function(object, correlation = FALSE, ...) {
  est <- object$coefficients
  se <- object$std.errs
  coe <- cbind(est, se, est / se, 2 * stats::pnorm(-abs(est / se)))
  dimnames(coe) <- list(names(est), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  ans <- list(call = object$call, terms = object$terms, residuals = object$residuals, coefficients = coe,
              aliased = is.na(est), sigma = object$sigma, df = c(object$rank, object$df.residual, length(est)),
              obj.function = object$obj.function, n.it = object$n.it, converged = object$converged,
              na.action = object$na.action)
  if (correlation) ans$correlation <- object$covar / outer(se[!is.na(se)], se[!is.na(se)])
  class(ans) <- "summary.Qlm"
  ans
}

#' @rdname summary.Qest
#' @export
print.summary.Qlm <- function(x, digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  resid <- x$residuals
  cat("Residuals:\n")
  if (length(resid) > 5L) {
    rq <- structure(zapsmall(stats::quantile(resid), digits + 1L), names = c("Min", "1Q", "Median", "3Q", "Max"))
    print(rq, digits = digits)
  } else print(resid, digits = digits)
  if (nsing <- x$df[3L] - x$df[1L]) cat("\nCoefficients: (", nsing, " not defined because of singularities)\n", sep = "")
  else cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  cat("\nsigma: ", format(signif(x$sigma, digits)), " (std. error of log(sigma): ",
      format(signif(attr(x$sigma, "std.err"), digits)), ")\n", sep = "")
  cat("Minimized loss function: ", format(x$obj.function, digits = digits + 3), "\n", sep = "")
  if (nzchar(mess <- stats::naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$correlation) && NCOL(x$correlation) > 1) {
    cat("\nCorrelation of Coefficients:\n")
    correl <- format(round(x$correlation, 2), nsmall = 2, digits = digits)
    correl[!lower.tri(correl)] <- ""
    print(correl[-1, -NCOL(correl), drop = FALSE], quote = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' Predictions from Qlm Models
#'
#' @param object an object of class `"Qlm"`.
#' @param newdata an optional data frame.
#' @param type `"response"` (the fitted mean, i.e. median, \eqn{x^\top\beta}) or
#'   `"quantile"` (\eqn{x^\top\beta + \sigma\Phi^{-1}(p)}).
#' @param se.fit logical: return the standard errors of the fitted mean?
#' @param p the quantile orders for `type = "quantile"`.
#' @param na.action how to treat missing values in `newdata`.
#' @param ... not used.
#' @return A vector, a list with items `fit` and `se.fit`, or (for
#'   `type = "quantile"`) a matrix with one column for each value of `p`.
#' @seealso [Qlm()]
#' @examples
#' set.seed(1)
#' x <- runif(100)
#' y <- 1 + 2 * x + rnorm(100)
#' fit <- Qlm(y ~ x)
#' predict(fit, newdata = data.frame(x = c(0, 0.5)), se.fit = TRUE)
#' predict(fit, newdata = data.frame(x = c(0, 0.5)), type = "quantile", p = c(0.1, 0.9))
#' @export
predict.Qlm <- function(object, newdata, type = c("response", "quantile"), se.fit = FALSE,
                        p = c(0.1, 0.25, 0.5, 0.75, 0.9), na.action = stats::na.pass, ...) {
  type <- match.arg(type)
  tt <- stats::delete.response(object$terms)
  X <- if (missing(newdata) || is.null(newdata)) stats::model.matrix(tt, object$model)
       else stats::model.matrix(tt, stats::model.frame(tt, newdata, na.action = na.action, xlev = object$xlevels),
                                contrasts.arg = object$contrasts)
  ok <- !is.na(object$coefficients)
  X <- X[, ok, drop = FALSE]
  fitv <- c(X %*% object$coefficients[ok])
  names(fitv) <- rownames(X)
  if (type == "quantile") {
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    out <- outer(fitv, object$sigma * stats::qnorm(p), "+")
    dimnames(out) <- list(rownames(X), paste0("p", p))
    return(out)
  }
  if (se.fit) list(fit = fitv, se.fit = sqrt(rowSums((X %*% object$covar) * X))) else fitv
}

#' @export
model.frame.Qlm <- function(formula, ...) formula$model

#' @export
model.matrix.Qlm <- function(object, ...) stats::model.matrix(object$terms, object$model, contrasts.arg = object$contrasts)

#' @export
nobs.Qlm <- function(object, ...) length(object$residuals)
