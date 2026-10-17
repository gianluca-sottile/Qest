# Methods for Qcoxph objects. All predictions use the fitted parametric model
# H(t | x) = H0(t) exp(x'beta), H0(t) = b(t)' exp(gamma).

#' @export
print.Qcoxph <- function(x, digits = max(1L, getOption("digits") - 3L), signif.stars = FALSE, ...) {
  if (!is.null(cl <- x$call)) { cat("Call:\n"); dput(cl); cat("\n") }
  coef <- x$coefficients
  se <- sqrt(diag(x$var))
  tmp <- cbind(coef, exp(coef), se, coef / se, stats::pchisq((coef / se)^2, 1, lower.tail = FALSE))
  dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "z", "p"))
  stats::printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE, signif.stars = signif.stars, ...)
  cat("\nObjective function = ", format(round(x$obj.function, 2)), " on ", sum(!is.na(coef)), " df\n", sep = "")
  cat("n = ", x$n, ", number of events = ", x$nevent, "\n", sep = "")
  if (length(x$na.action)) cat("   (", stats::naprint(x$na.action), ")\n", sep = "")
  invisible(x)
}

#' Summarizing Qcoxph Fits
#'
#' @param object an object of class `"Qcoxph"`.
#' @param conf.int the level of the confidence intervals for the hazard ratios.
#' @param scale a vector of scale factors for the coefficients.
#' @param x an object of class `"summary.Qcoxph"`.
#' @param digits,signif.stars printing options.
#' @param ... not used.
#' @return An object of class `"summary.Qcoxph"`, with the table of
#'   coefficients (`coefficients`), the confidence intervals of the hazard ratios
#'   (`conf.int`), the concordance index and its standard error (`concordance`),
#'   and the loss (`obj.function`).
#' @seealso [Qcoxph()]
#' @export
summary.Qcoxph <- function(object, conf.int = 0.95, scale = 1, ...) {
  beta <- object$coefficients * scale
  se <- sqrt(diag(object$var)) * scale
  rval <- list(call = object$call, na.action = object$na.action, n = object$n, nevent = object$nevent,
               obj.function = object$obj.function, type = object$internal$type)
  tmp <- cbind(beta, exp(beta), se, beta / se, stats::pchisq((beta / se)^2, 1, lower.tail = FALSE))
  dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
  rval$coefficients <- tmp
  if (conf.int) {
    z <- stats::qnorm((1 + conf.int) / 2)
    tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se), exp(beta + z * se))
    dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)", paste0("lower .", round(100 * conf.int, 2)),
                                          paste0("upper .", round(100 * conf.int, 2))))
    rval$conf.int <- tmp
  }
  if (!is.null(object$concordance)) {
    rval$concordance <- object$concordance[c("concordance", "std")]
    names(rval$concordance) <- c("C", "se(C)")
  }
  class(rval) <- "summary.Qcoxph"
  rval
}

#' @rdname summary.Qcoxph
#' @export
print.summary.Qcoxph <- function(x, digits = max(getOption("digits") - 3, 3),
                                 signif.stars = getOption("show.signif.stars"), ...) {
  if (!is.null(x$call)) { cat("Call:\n"); dput(x$call); cat("\n") }
  cat("  n = ", x$n, ", number of events = ", x$nevent, "\n", sep = "")
  if (length(x$na.action)) cat("   (", stats::naprint(x$na.action), ")\n", sep = "")
  cat("\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  if (!is.null(x$conf.int)) { cat("\n"); print(x$conf.int, digits = digits) }
  cat("\n")
  if (!is.null(x$concordance))
    cat("Concordance = ", format(round(x$concordance[1], 3)), " (se = ", format(round(x$concordance[2], 3)), ")\n",
        sep = "")
  cat("Objective function = ", format(round(x$obj.function, 2), digits = digits),
      if (x$type != "u") " (not the function being minimized)", "\n", sep = "")
  invisible(x)
}

#' @export
vcov.Qcoxph <- function(object, ...) object$var

#' @export
confint.Qcoxph <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients
  se <- sqrt(diag(object$var))
  ci_table(cf, se, parm, level)
}

#' @export
nobs.Qcoxph <- function(object, ...) object$n

#' @export
model.matrix.Qcoxph <- function(object, ...) object$x

# Design matrix (without intercept) for new data.
cox_newx <- function(object, newdata) {
  if (missing(newdata) || is.null(newdata)) return(object$x)
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = object$xlevels)
  X <- stats::model.matrix(tt, mf)
  X[, colnames(object$x), drop = FALSE]
}

# H0(t) and its gradient w.r.t. (beta, gamma), for the observed covariates.
cox_H <- function(object, times, X) {
  ok <- object$internal$ok
  beta <- object$coefficients[ok]
  gamma <- object$internal$gamma
  b <- plfcox(times, object$internal$knots)
  H0 <- c(b %*% exp(gamma))
  eta <- c(X[, ok, drop = FALSE] %*% beta)
  # H[t, i] = H0(t) exp(eta_i); gradient: dH/dbeta = H x_i, dH/dgamma_l = exp(eta_i) b_l(t) exp(gamma_l)
  list(H0 = H0, eta = eta, b = b, egamma = exp(gamma), X = X[, ok, drop = FALSE])
}

# Variance of H(t | x_i) by the delta method, for all t and one x.
cox_varH <- function(cH, i) {
  covar <- cH$covar
  H <- cH$H0 * exp(cH$eta[i])
  G <- cbind(H %o% cH$X[i, ], exp(cH$eta[i]) * t(t(cH$b) * cH$egamma))
  rowSums((G %*% covar) * G)
}

#' Baseline Cumulative Hazard of a Qcoxph Model
#'
#' @param fit an object of class `"Qcoxph"`.
#' @param centered logical: if `TRUE`, the cumulative hazard at the mean of the
#'   covariates; otherwise at \eqn{x = 0}.
#' @param se.fit logical: compute the variance of the estimates (delta method)?
#' @param times the times at which \eqn{H_0} is computed (by default, the
#'   observed times).
#' @return A data frame with columns `hazard` (the cumulative hazard), `time`
#'   and, if `se.fit = TRUE`, `varhaz`.
#' @details The baseline cumulative hazard of the fitted model is piecewise
#'   linear, \eqn{H_0(t) = \sum_l \exp(\gamma_l) b_l(t)}. Unlike
#'   [survival::basehaz()], this is not a step function.
#' @seealso [Qcoxph()], [survfit.Qcoxph()]
#' @examples
#' library(survival)
#' set.seed(1)
#' x <- runif(200)
#' t <- rexp(200, exp(x))
#' cens <- runif(200, 0, 3)
#' fit <- Qcoxph(Surv(pmin(t, cens), t <= cens) ~ x)
#' head(basehaz.Qcoxph(fit, se.fit = TRUE))
#' @export
basehaz.Qcoxph <- function(fit, centered = TRUE, se.fit = FALSE, times) {
  if (!inherits(fit, "Qcoxph")) stop("'fit' must be an object of class \"Qcoxph\"", call. = FALSE)
  if (missing(times)) times <- sort(unique(fit$internal$y))
  X <- matrix(if (centered) fit$means else 0, 1, length(fit$means), dimnames = list(NULL, names(fit$means)))
  cH <- cox_H(fit, times, X)
  cH$covar <- fit$internal$covar
  out <- data.frame(hazard = cH$H0 * exp(cH$eta[1]), time = times)
  if (se.fit) out$varhaz <- cox_varH(cH, 1)
  out
}

#' Predictions from a Qcoxph Model
#'
#' @param object an object of class `"Qcoxph"`.
#' @param newdata an optional data frame. For `type = "expected"` and
#'   `type = "survival"` it must contain the response variables.
#' @param type the type of prediction: the linear predictor (`"lp"`, centered at
#'   the mean of the covariates, as in [survival::coxph()]), the relative risk
#'   (`"risk"`, \eqn{\exp(lp)}), the expected number of events
#'   \eqn{H(y \mid x)} (`"expected"`, minus \eqn{H(z \mid x)} with truncation), the
#'   survival probability \eqn{\exp(-expected)} (`"survival"`), or the quantiles of
#'   the time-to-event (`"quantile"`).
#' @param se.fit logical: compute standard errors (delta method)? Not available
#'   for `type = "quantile"`.
#' @param p the quantile orders for `type = "quantile"`.
#' @param ... not used.
#' @return A vector, or a list with items `fit` and `se.fit`; for
#'   `type = "quantile"`, a matrix with one column for each value of `p`.
#' @seealso [Qcoxph()], [survfit.Qcoxph()]
#' @examples
#' library(survival)
#' set.seed(1)
#' dat <- data.frame(x = runif(200))
#' t <- rexp(200, exp(dat$x)); cens <- runif(200, 0, 3)
#' dat$y <- pmin(t, cens); dat$d <- as.numeric(t <= cens)
#' fit <- Qcoxph(Surv(y, d) ~ x, data = dat)
#' predict(fit, newdata = data.frame(x = c(0, 0.5, 1)), type = "risk", se.fit = TRUE)
#' predict(fit, newdata = data.frame(x = c(0, 1)), type = "quantile", p = c(0.25, 0.5))
#' @export
predict.Qcoxph <- function(object, newdata, type = c("lp", "risk", "expected", "survival", "quantile"),
                           se.fit = FALSE, p = c(0.25, 0.5, 0.75), ...) {
  type <- match.arg(type)
  ok <- object$internal$ok
  X <- cox_newx(object, if (missing(newdata)) NULL else newdata)
  beta <- object$coefficients[ok]
  Xc <- sweep(X[, ok, drop = FALSE], 2, object$means[ok])
  V <- object$var[ok, ok, drop = FALSE]
  if (type %in% c("lp", "risk")) {
    lp <- c(Xc %*% beta)
    se <- sqrt(rowSums((Xc %*% V) * Xc))
    names(lp) <- names(se) <- rownames(X)
    if (type == "risk") { fitv <- exp(lp); se <- fitv * se } else fitv <- lp
    return(if (se.fit) list(fit = fitv, se.fit = se) else fitv)
  }
  if (type == "quantile") {
    if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
    eta <- c(X[, ok, drop = FALSE] %*% beta)
    knots <- object$internal$knots
    H0k <- c(plfcox(knots, knots) %*% exp(object$internal$gamma))
    out <- sapply(p, function(pp) {
      h <- -log(1 - pp) * exp(-eta)
      q <- stats::approx(H0k, knots, xout = h)$y
      K <- length(knots); up <- h > H0k[K]   # H0 is linear after the last knot
      q[up] <- knots[K] + (h[up] - H0k[K]) / exp(object$internal$gamma[K - 1])
      q
    })
    out <- matrix(out, nrow(X))
    dimnames(out) <- list(rownames(X), paste0("p", p))
    return(out)
  }
  # expected / survival
  if (missing(newdata)) { y <- object$internal$y; z <- object$internal$z; type0 <- object$internal$type }
  else {
    mf <- stats::model.frame(object$terms, newdata, na.action = stats::na.pass, xlev = object$xlevels)
    r <- get_response(stats::model.response(mf))
    y <- r$y; z <- r$z; type0 <- r$type
  }
  cH <- cox_H(object, y, X)
  cH$covar <- object$internal$covar
  b <- plfcox(y, object$internal$knots)
  bz <- if (type0 == "ct") plfcox(pmax(z, object$internal$knots[1]), object$internal$knots) else 0 * b
  eg <- exp(object$internal$gamma)
  e <- exp(cH$eta)
  Hy <- c((b - bz) %*% eg) * e
  G <- cbind(Hy * cH$X, e * t(t(b - bz) * eg))
  se <- sqrt(rowSums((G %*% cH$covar) * G))
  if (type == "survival") { fitv <- exp(-Hy); se <- fitv * se } else fitv <- Hy
  names(fitv) <- names(se) <- rownames(X)
  if (se.fit) list(fit = fitv, se.fit = se) else fitv
}

#' Residuals of a Qcoxph Model
#'
#' @param object an object of class `"Qcoxph"`.
#' @param type `"martingale"` (\eqn{d_i - H(y_i \mid x_i)}) or `"deviance"`.
#' @param ... not used.
#' @return A vector of residuals.
#' @seealso [Qcoxph()]
#' @export
residuals.Qcoxph <- function(object, type = c("martingale", "deviance"), ...) {
  type <- match.arg(type)
  rr <- object$residuals
  if (type == "deviance") {
    d <- object$internal$d
    rr <- sign(rr) * sqrt(-2 * (rr + ifelse(d == 0, 0, d * log(d - rr))))
  }
  rr
}

#' Survival Curves from a Qcoxph Model
#'
#' Computes the survival curves \eqn{S(t \mid x) = \exp\{-H_0(t)\exp(x^\top\beta)\}}
#' of the fitted model, with pointwise confidence intervals based on the
#' delta method.
#'
#' @param formula an object of class `"Qcoxph"`.
#' @param newdata a data frame with the covariates of the curves. By default, one
#'   curve at the mean of the covariates.
#' @param se.fit logical: compute the standard errors?
#' @param conf.int the level of the confidence intervals.
#' @param conf.type the transformation used for the confidence intervals:
#'   `"log"` (default, symmetric on the scale of \eqn{H}), `"log-log"`,
#'   `"plain"` or `"none"`.
#' @param times the times (by default, the observed times).
#' @param ... not used.
#' @return An object of class `"survfit"` (see [survival::survfit.object]): the
#'   methods `print`, `plot`, `lines` and `summary` of 'survival' can be used.
#' @seealso [Qcoxph()], [predict.Qcoxph()]
#' @examples
#' library(survival)
#' set.seed(1)
#' dat <- data.frame(x = runif(200))
#' t <- rexp(200, exp(dat$x)); cens <- runif(200, 0, 3)
#' dat$y <- pmin(t, cens); dat$d <- as.numeric(t <= cens)
#' fit <- Qcoxph(Surv(y, d) ~ x, data = dat)
#' sf <- survfit(fit, newdata = data.frame(x = c(0, 1)))
#' plot(sf, col = 1:2, conf.int = TRUE)
#' @importFrom survival survfit
#' @export
survfit.Qcoxph <- function(formula, newdata, se.fit = TRUE, conf.int = 0.95,
                           conf.type = c("log", "log-log", "plain", "none"), times, ...) {
  object <- formula
  conf.type <- match.arg(conf.type)
  if (!se.fit) conf.type <- "none"
  X <- if (missing(newdata)) matrix(object$means, 1, dimnames = list(NULL, names(object$means)))
       else cox_newx(object, newdata)
  y <- object$internal$y; d <- object$internal$d; w <- object$internal$w
  if (missing(times)) times <- sort(unique(y))
  times <- sort(times)
  cH <- cox_H(object, times, X)
  cH$covar <- object$internal$covar
  H <- outer(cH$H0, exp(cH$eta))
  S <- exp(-H)
  seH <- vapply(seq_len(nrow(X)), function(i) sqrt(pmax(cox_varH(cH, i), 0)), numeric(length(times)))
  seH <- matrix(seH, length(times))
  z <- stats::qnorm((1 + conf.int) / 2)
  lower <- upper <- NULL
  if (conf.type == "log") { lower <- exp(-(H + z * seH)); upper <- pmin(exp(-(H - z * seH)), 1) }
  else if (conf.type == "log-log") {
    s <- seH / pmax(H, 1e-300)
    lower <- exp(-H * exp(z * s)); upper <- exp(-H * exp(-z * s))
  } else if (conf.type == "plain") { lower <- pmax(S - z * S * seH, 0); upper <- pmin(S + z * S * seH, 1) }
  nrisk <- vapply(times, function(t) sum(w[y >= t]), 1)
  nevent <- vapply(times, function(t) sum(w[y == t & d == 1]), 1)
  ncens <- vapply(times, function(t) sum(w[y == t & d == 0]), 1)
  out <- list(n = object$n, time = times, n.risk = nrisk, n.event = nevent, n.censor = ncens,
              surv = drop_one(S), cumhaz = drop_one(H), std.err = if (se.fit) drop_one(seH),
              lower = drop_one(lower), upper = drop_one(upper),
              conf.type = conf.type, conf.int = if (conf.type != "none") conf.int, logse = TRUE,
              type = "right", call = match.call())
  out <- out[!vapply(out, is.null, TRUE)]
  class(out) <- "survfit"
  out
}

drop_one <- function(M) if (is.null(M)) NULL else if (ncol(M) == 1) c(M) else M
