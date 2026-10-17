#' Quantile-Based Estimators (Q-Estimators)
#'
#' Fits a parametric model by minimizing the integrated loss of quantile
#' regression (Sottile and Frumento, 2022). The model is specified through its
#' quantile function \eqn{Q(\tau \mid x, \theta)}, either written by the user
#' or chosen among the built-in families (see [Qfamily]). Censored and
#' truncated data are allowed.
#'
#' @param formula a two-sided formula of the form `y ~ x`. The response can be
#'   a numeric vector, a right-censored `survival::Surv(time, event)` object, or
#'   a left-truncated and right-censored `survival::Surv(start, stop, event)`
#'   object. With a user-defined `Q`, the right-hand side of the formula is only
#'   used to define the covariates of the starting points: write `y ~ 1` if no
#'   covariates are used by `Q`.
#' @param Q a quantile function `Q(theta, tau, data)` (see "Details"), or a
#'   [Qfamily] (`Qnorm`, `Qgamma`, `Qpois`, `Qunif`).
#' @param weights an optional vector of non-negative weights. They are treated as
#'   sampling weights: the (sandwich) standard errors do not depend on their scale.
#' @param start the starting points. With a user-defined `Q`, `start` is
#'   required, because its length defines the number of parameters: use `NA`
#'   for the parameters without a good starting value (with any `NA`, the
#'   starting points are estimated, see "Details"). With a `Qfamily`, `start`
#'   is optional.
#' @param data an optional data frame containing the variables in the model. It
#'   is also the argument `data` of `Q`.
#' @param subset,na.action as in [stats::lm()].
#' @param ntau the number of points of the grid of quantile orders used for
#'   numerical integration (default 199).
#' @param wtau an optional weighting function \eqn{w(\tau)}, e.g. [wtrunc()].
#'   Its first argument must be `tau`; other arguments are passed through `...`.
#' @param control a list of operational parameters created by [Qest.control()].
#' @param ... additional named arguments of `wtau` or `Q`.
#'
#' @details
#' The Q-estimator minimizes
#' \deqn{L(\theta) = \sum_i \int_0^1 w(\tau)\rho_\tau(y_i - Q(\tau \mid x_i, \theta))\,d\tau,}
#' where \eqn{\rho_\tau(u) = u(\tau - I(u < 0))}. It is consistent and
#' asymptotically normal, and its robustness can be tuned by choosing
#' \eqn{w(\tau)} (Sottile and Frumento, 2022). For censored and truncated data,
#' the estimating equations of Appendix A of the paper are solved.
#'
#' **User-defined quantile functions.** `Q(theta, tau, data)` must have these
#' three arguments, in this order (other arguments can follow and are passed
#' through `...`), and must return, for each row of `data`,
#' \eqn{Q(\tau \mid x_i, \theta)}. The argument `tau` is either a matrix with
#' one row for each row of `data` (then `Q` returns a matrix of the same size),
#' or a vector with one value for each row of `data` (then `Q` returns a
#' vector). `data` is the data frame used to fit the model (after removing
#' missing values, `subset` and observations with zero weight). `Q` must be
#' increasing in `tau` and twice differentiable in `theta`; it is evaluated on a grid of `ntau` quantile orders in
#' \eqn{[5 \cdot 10^{-4}, 1 - 5 \cdot 10^{-4}]}, and its derivatives with respect
#' to `theta` are computed by central differences.
#'
#' **Starting points.** If `start` contains `NA` (or `control$restart = TRUE`),
#' a nonparametric estimate \eqn{\hat F(y \mid x)} is obtained from a
#' piecewise-constant hazard model that uses the covariates in `formula`, and the
#' starting points minimize \eqn{\sum_i (y_i - Q(\hat F(y_i \mid x_i) \mid x_i, \theta))^2}
#' (Sottile and Frumento, 2022, Algorithm 1), searched from the values in `start`
#' (with `NA` replaced by zero): with poorly scaled data, rough starting values for
#' all parameters are recommended.
#'
#' **Algorithm.** The estimating equations are solved by Newton-Raphson with
#' backtracking, switching to Levenberg-Marquardt steps when the Newton step
#' fails. The Jacobian is the exact derivative of the discretized estimating
#' equations. Standard errors are computed with the sandwich formula.
#'
#' @return An object of class `"Qest"`, a list with items
#' \item{coefficients}{the estimated parameters.}
#' \item{std.errs}{their standard errors.}
#' \item{covar}{the estimated covariance matrix.}
#' \item{obj.function}{the value of the loss at the solution (for censored or
#'   truncated data, the loss of Appendix A, which is not the function whose
#'   gradient is solved).}
#' \item{ee, jacobian}{the estimating equations (divided by the sample size)
#'   and their Jacobian at the solution.}
#' \item{CDF, PDF}{the fitted \eqn{F(y_i \mid x_i)} and \eqn{f(y_i \mid x_i)}.}
#' \item{converged, n.it}{convergence of the algorithm and number of iterations.}
#' \item{internal}{quantities used by the methods.}
#' and the usual `call`, `formula`, `terms`, `data`, `na.action`, `xlevels`,
#' `contrasts`, `control`.
#'
#' @references Sottile, G. and Frumento, P. (2022). Robust estimation and
#'   regression with parametric quantile functions. *Computational Statistics
#'   and Data Analysis*, 171, 107471. \doi{10.1016/j.csda.2022.107471}
#' @seealso [Qfamily], [Qlm()], [Qcoxph()], [summary.Qest()], [wtrunc()],
#'   [invQ()].
#' @examples
#' # Normal linear model
#' Qlinmod <- function(theta, tau, data) {
#'   X <- model.matrix(~ x1 + x2, data = data)
#'   qnorm(tau, c(X %*% theta[-1]), exp(theta[1]))
#' }
#' set.seed(1)
#' n <- 100
#' dat <- data.frame(x1 = rnorm(n), x2 = runif(n, 0, 3))
#' dat$y <- 4 + dat$x1 + 2 * dat$x2 + rnorm(n)
#' m1 <- Qest(y ~ x1 + x2, Q = Qlinmod, start = rep(NA, 4), data = dat)
#' m2 <- Qest(y ~ x1 + x2, Q = Qnorm, data = dat)    # the same model, as a Qfamily
#' m3 <- Qlm(y ~ x1 + x2, data = dat)                # ... or with Qlm (fastest)
#' summary(m1)
#' cbind(user = coef(m1), family = coef(m2))
#'
#' \donttest{
#' # A Weibull proportional hazards model with censored data
#' QWeibPH <- function(theta, tau, data) {
#'   shape <- exp(theta[1])
#'   X <- model.matrix(~ x, data = data)
#'   qweibull(tau, shape = shape, scale = exp(-c(X %*% theta[-1]) / shape))
#' }
#' dat <- data.frame(x = runif(n, 0, 3))
#' t <- QWeibPH(c(log(2), -0.5, 1), runif(n), dat)
#' cens <- runif(n, 0, 2)
#' dat$y <- pmin(t, cens); dat$d <- as.numeric(t <= cens)
#' m4 <- Qest(survival::Surv(y, d) ~ x, Q = QWeibPH, start = rep(NA, 3), data = dat)
#' summary(m4)
#'
#' # Gamma model with more weight around the median
#' y <- rgamma(n, shape = 3, scale = 2)
#' m5 <- Qest(y ~ 1, Q = Qgamma, wtau = function(tau, h) dnorm((tau - 0.5) / h), h = 0.2)
#' coef(m5)
#' }
#' @export
Qest <- function(formula, Q, weights, start, data, subset, na.action, ntau = 199, wtau = NULL,
                 control = Qest.control(), ...) {
  cl <- match.call()
  Q <- resolve_Q(Q, parent.frame())
  family <- inherits(Q, "Qfamily")

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mf0 <- mf

  w <- check_weights(stats::model.weights(mf), nrow(mf))
  if (any(w == 0)) {
    warning("observations with zero weight are dropped", call. = FALSE)
    mf <- mf[w > 0, , drop = FALSE]
    w <- w[w > 0]
  }
  w <- w / mean(w)
  if ((nobs <- nrow(mf)) == 0) stop("zero non-NA cases", call. = FALSE)

  # the data passed to Q: the rows of 'data' used in the model frame
  if (missing(data)) Qdata <- mf
  else {
    if (!is.data.frame(data)) data <- as.data.frame(data)
    Qdata <- data[match(rownames(mf), rownames(data)), , drop = FALSE]
  }

  resp <- get_response(stats::model.response(mf))
  y <- resp$y; z <- resp$z; d <- resp$d; type <- resp$type
  x <- if (!stats::is.empty.model(mt)) stats::model.matrix(mt, mf) else matrix(1, nobs, 1, dimnames = list(NULL, "(Intercept)"))

  dots <- split_dots(list(...), wtau, if (!family) Q)
  if (!is.null(wtau) && !is.function(wtau)) stop("'wtau' must be a function", call. = FALSE)
  if (!is.null(wtau) && names(formals(wtau))[1] != "tau") stop("the first argument of 'wtau' must be 'tau'", call. = FALSE)
  control <- do.call(Qest.control, control)

  fit <- if (family) {
    off <- stats::model.offset(mf)
    if (!is.null(off) && Q$family$family != "poisson") stop("offsets are only allowed with Qpois", call. = FALSE)
    qest_family_fit(Q, z, y, d, x, w, type, if (missing(start)) NULL else start, ntau, wtau, dots$wtau,
                    control, offset = off)
  } else {
    if (missing(start) || is.null(start))
      stop("please provide 'start': its length defines the number of parameters.\n",
           "Use NA for unknown values, e.g., start = c(NA, NA, NA).", call. = FALSE)
    qest_generic_fit(Q, dots$Q, Qdata, z, y, d, x, w, type, start, ntau, wtau, dots$wtau, control)
  }

  CDF <- fit$CDF; PDF <- fit$PDF
  names(CDF) <- names(PDF) <- rownames(mf)
  out <- list(coefficients = fit$theta, std.errs = fit$se, covar = fit$covar,
              obj.function = fit$loss, ee = fit$ee, jacobian = fit$jacobian,
              CDF = CDF, PDF = PDF, converged = fit$converged, n.it = fit$n.it,
              internal = list(z = z, y = y, d = d, w = w, Q = fit$Q, type = type, tau = fit$tau,
                              family = fit$family, Qoptions = dots$Q, ntau = ntau, wtau = wtau,
                              wtauoptions = dots$wtau, x = x, offset = fit$offset),
              formula = formula, terms = mt, data = Qdata, control = control,
              na.action = attr(mf0, "na.action"), xlevels = stats::.getXlevels(mt, mf),
              contrasts = attr(x, "contrasts"), call = cl)
  class(out) <- "Qest"
  if (!fit$converged) warning("the algorithm did not converge", call. = FALSE)
  out
}

#' Control Parameters of Qest, Qlm and Qcoxph
#'
#' @param tol convergence tolerance on the norm of the estimating equations
#'   (divided by the sample size).
#' @param maxit maximum number of iterations (default `50 + 25 * npar`).
#' @param display logical: print the iterations?
#' @param restart logical: if `TRUE`, the starting points of [Qest()] with a
#'   user-defined `Q` are re-estimated even if `start` has no missing values.
#' @param deriv.step relative step of the central differences used for the
#'   derivatives of a user-defined `Q` with respect to `theta`: the step for
#'   \eqn{\theta_j} is `deriv.step * max(1, abs(theta_j))`.
#' @param ... arguments of 'Qest' 1.0.x that were removed (`safeit`, `alpha0`):
#'   they give an error.
#' @return A list with the control parameters.
#' @seealso [Qest()], [Qlm()], [Qcoxph()]
#' @examples
#' Qest.control(tol = 1e-10, display = TRUE)
#' @export
Qest.control <- function(tol = 1e-8, maxit = NULL, display = FALSE, restart = FALSE,
                         deriv.step = 1e-4, ...) {
  check_removed(list(...), "Qest.control", c(safeit = "", alpha0 = ""))
  if (!is.numeric(tol) || tol <= 0) stop("'tol' must be positive", call. = FALSE)
  if (!is.null(maxit) && (!is.numeric(maxit) || maxit < 1)) stop("'maxit' must be a positive number", call. = FALSE)
  if (!is.numeric(deriv.step) || deriv.step <= 0 || deriv.step > 0.1)
    stop("'deriv.step' must be in (0, 0.1]", call. = FALSE)
  list(tol = tol, maxit = maxit, display = display, restart = restart, deriv.step = deriv.step)
}

#' @rdname Qest.control
#' @export
Qcoxph.control <- function(tol = 1e-8, maxit = NULL, display = FALSE, ...) {
  check_removed(list(...), "Qcoxph.control", c(safeit = "", alpha0 = ""))
  Qest.control(tol = tol, maxit = maxit, display = display)
}

# Error for arguments of version 1.0.x that were removed.
check_removed <- function(dots, fun, removed = character(0)) {
  if (!length(dots)) return(invisible())
  nm <- names(dots)
  if (is.null(nm) || any(nm == "")) stop("unused unnamed argument(s) in ", fun, "()", call. = FALSE)
  rem <- intersect(nm, names(removed))
  if (length(rem))
    stop("argument(s) ", paste(sQuote(rem), collapse = ", "), " of ", fun,
         "() were removed in Qest 2.0.0 (the algorithm no longer uses them)", call. = FALSE)
  stop("unused argument(s) in ", fun, "(): ", paste(nm, collapse = ", "), call. = FALSE)
}

check_weights <- function(w, n) {
  if (is.null(w)) return(rep.int(1, n))
  if (!is.numeric(w)) stop("'weights' must be a numeric vector", call. = FALSE)
  if (any(!is.finite(w)) || any(w < 0)) stop("'weights' must be non-negative and finite", call. = FALSE)
  if (all(w == 0)) stop("all weights are zero", call. = FALSE)
  w
}

# Response: numeric, Surv(time, event) or Surv(start, stop, event).
get_response <- function(y) {
  n <- NROW(y)
  type <- if (inherits(y, "Surv")) attr(y, "type") else NULL
  if (is.null(type)) {
    if (!is.numeric(y)) stop("the response must be numeric or a 'Surv' object", call. = FALSE)
    y <- as.numeric(y)
    return(list(y = y, z = rep.int(-Inf, n), d = rep.int(1, n), type = "u"))
  }
  y <- unclass(y)
  if (type == "right") {
    out <- list(y = y[, 1], z = rep.int(-Inf, n), d = y[, 2])
    out$type <- if (any(out$d == 0)) "c" else "u"
  } else if (type == "counting") {
    out <- list(z = y[, 1], y = y[, 2], d = y[, 3])
    if (any(out$z >= out$y)) stop("the truncation time must be smaller than the observed time", call. = FALSE)
    out$type <- "ct"
    # (truncation matters even if all z < min(y), because F(z | x) > 0 in general)
    if (all(!is.finite(out$z))) {
      out$type <- if (any(out$d == 0)) "c" else "u"
      out$z <- rep.int(-Inf, n)
    }
  } else stop("only right-censored and left-truncated data are supported", call. = FALSE)
  if (!any(out$d == 1)) stop("all data are censored", call. = FALSE)
  out
}

resolve_Q <- function(Q, env) {
  if (missing(Q)) stop("'Q' is missing", call. = FALSE)
  if (is.character(Q)) {
    if (!Q %in% c("Qnorm", "Qgamma", "Qpois", "Qunif")) stop("'Qfamily' not recognized", call. = FALSE)
    Q <- match.fun(Q)
  }
  if (inherits(Q, "Qfamily")) return(Q)
  if (!is.function(Q)) stop("'Q' must be a function or a 'Qfamily'", call. = FALSE)
  for (f in list(Qnorm, Qgamma, Qpois, Qunif)) if (identical(Q, f)) return(Q())
  argQ <- names(formals(Q))
  if (length(argQ) < 3 || any(argQ[1:3] != c("theta", "tau", "data")))
    stop("Q must have arguments 'theta', 'tau', 'data', in this order, followed by other optional arguments",
         call. = FALSE)
  Q
}

# Fit with a user-defined Q.
qest_generic_fit <- function(Q, Qoptions, data, z, y, d, x, w, type, start, ntau, wtau, wtauoptions, control) {
  n <- length(y)
  attr(Q, "Qoptions") <- Qoptions
  npar <- length(start)
  tau <- buildTau(ntau, wtau, n, wtauoptions)
  test <- suppressWarnings(try(callQ(Q, ifelse(is.na(start), 0, start), tau$TAU[, 1:9, drop = FALSE], data),
                               silent = TRUE))
  if (inherits(test, "try-error")) stop("Q(theta, tau, data) returned an error: ", attr(test, "condition")$message,
                                        call. = FALSE)
  if (!is.matrix(test) || any(dim(test) != c(n, 9)))
    stop("Q(theta, tau, data) must return a matrix with nrow(data) rows when tau is a matrix", call. = FALSE)

  theta0 <- start_generic(z, y, d, x, w, Q, data, start, control)
  ctx <- list(Q = Q, data = data, y = y, z = z, d = d, w = w, type = type, tau = tau, h = control$deriv.step)
  eefun <- function(th, J) generic_ee(th, ctx, J)
  if (is.null(eefun(theta0, FALSE)))
    stop("Q(start, tau, data) does not return finite, increasing values: please provide better starting points",
         call. = FALSE)
  fit <- qest_solve(theta0, eefun, n, control)
  G <- fit$G
  theta <- fit$theta
  # The analytical Jacobian assumes that Q is twice differentiable in theta: compare
  # it with finite differences of the estimating equations.
  hJ <- 1e-3 * pmax(abs(theta), 1)
  gp <- lapply(seq_len(npar), function(j) {
    e <- replace(numeric(npar), j, hJ[j])
    a <- eefun(theta + e, FALSE); b <- eefun(theta - e, FALSE)
    if (is.null(a) || is.null(b)) NULL else (a$g - b$g) / (2 * hJ[j])
  })
  if (!any(vapply(gp, is.null, TRUE))) {
    Jfd <- do.call(cbind, gp)
    if (max(abs(G$J - Jfd)) > 0.05 * max(abs(Jfd))) {
      warning("the Jacobian computed from the derivatives of Q differs from a numerical Jacobian: ",
              "is Q twice differentiable with respect to theta? The numerical Jacobian is used ",
              "to compute the standard errors", call. = FALSE)
      G$J <- Jfd
    }
  }
  nms <- names(start)
  if (is.null(nms) || any(nms == "")) nms <- paste0("Theta", seq_len(npar))
  covar <- tryCatch(sandwich(G$J, G$g.i, w), error = function(e) NULL)
  if (is.null(covar)) {
    warning("the Jacobian is singular: is the model identifiable?", call. = FALSE)
    covar <- matrix(NA_real_, npar, npar)
  }
  names(theta) <- nms
  dimnames(covar) <- list(nms, nms)
  J <- G$J / n; dimnames(J) <- list(nms, nms)
  ee <- G$g / n; names(ee) <- nms
  loss <- generic_loss(G$Q0, y, z, d, w, type, tau, G$Fy$F, if (type == "ct") G$Fz$F)
  attr(loss, "npar") <- npar
  if (type != "u") attr(loss, "message") <- "Not the function being minimized"
  list(theta = theta, se = sqrt(pmax(diag(covar), 0)), covar = covar, loss = loss, ee = ee, jacobian = J,
       CDF = G$Fy$F, PDF = G$fy, converged = fit$converged, n.it = fit$n.it, Q = Q, tau = tau, family = NULL)
}
