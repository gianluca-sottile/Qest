# Qest 2.0.0

A major revision. Most estimates agree with version 1.0.2 up to small numerical
differences, but several inaccuracies were corrected (see "Numerical accuracy")
and the API was cleaned.

## Numerical accuracy

* User-defined quantile functions: the derivatives of `Q` with respect to `theta`
  are now computed by central differences with a small step (they were forward
  differences with a large, data-driven step), the integrals on the tails of the
  grid of quantile orders are no longer ignored, and the Jacobian is the exact
  derivative of the discretized estimating equations. The estimates now agree with
  the exact (family) estimators within about 0.001 standard errors, and the
  standard errors within 0.5% (in 1.0.2: about 1% of a standard error, and up to
  5-10% in the standard errors).
* With censored or truncated data, the Jacobian of the estimating equations (which
  is not symmetric) was used transposed in `Qest()` (user-defined `Q`, `Qnorm`,
  `Qgamma`) and `Qcoxph()`: this gave wrong standard errors (e.g., 20% too small for
  `log(sigma)` in a truncated normal model) and slowed down convergence. Fixed.
* `Qgamma`: the integral on the upper tail of the grid was ignored, which biased the
  estimates by about 1.5% of their standard error. Fixed.
* `Qpois`: F(y | x) was computed exactly while the integrals were interpolated on the
  approximated quantiles, so that the estimating equations were not smooth (the
  algorithm stopped up to 0.1 standard errors away from the solution). The
  estimating equations are now consistent with the approximated quantile function
  and converge quadratically; the estimates agree with an exact continuous Poisson
  quantile function within 0.001 standard errors.
* `Qnorm` with left-truncated data: the truncation times were not rescaled.
* `Qunif` with `wtau`: the last step of the integration grid was zero.

## Algorithm

* All models use the same solver: Newton-Raphson with backtracking on the norm of
  the estimating equations, switching to Levenberg-Marquardt steps when needed, with
  a limit on the step length. The stochastic and gradient-search preliminary phases
  were removed: results no longer depend on the random seed (in 1.0.2, `Qlm()` and
  `Qest()` used random sub-samples).
* Starting points: a piecewise-constant hazard model written in the package
  replaces the internal functions of 'pch'; for `Qcoxph()`, the partial likelihood
  and the Breslow estimator are computed in the package (no internal functions of
  'survival'), and the knot-moving step copied from 'segmented' was removed (it was
  never applied to the default knots, and applied inconsistently to user knots).
* The Fortran code (Barrodale-Roberts algorithm, used for starting points) was
  replaced by an R implementation of the Frisch-Newton algorithm: the package no
  longer needs compilation.

## API changes

* `Qest()` gains `subset` and `na.action`; the rows of `data` passed to `Q` are now
  those used in the model (in 1.0.2, rows with missing values or zero weights were
  passed to `Q`). Unknown arguments in `...` are an error (they were ignored).
* `Qest.control()` and `Qcoxph.control()`: the arguments `safeit` and `alpha0` were
  removed (they give an error); new argument `deriv.step`.
* `Qlm()` returns an object of class `c("Qlm", "Qest")` (no longer `"lm"`), with new
  items `sigma` and `covar.full`; `summary.Qlm()` uses normal (asymptotic) tests and
  no longer reports least squares R-squared and F statistics; new `predict.Qlm()`.
  `Qlm(y ~ 0)` gives an informative error.
* `Qcoxph()` returns an object of class `c("Qcoxph", "Qest")` (no longer `"coxph"`,
  whose methods did not apply); gains `subset` and `na.action`. `survfit()`,
  `predict()`, `residuals()` and `basehaz.Qcoxph()` were rewritten from the fitted
  parametric model (they copied about 1000 lines of 'survival', including a call to
  a C routine of 'survival'); `survfit()` returns a standard `"survfit"` object and
  `predict()` gains `type = "quantile"`. Martingale residuals account for left
  truncation.
* New methods: `predict.Qest()` (quantiles and CDF), `plot.Qest()` (uniform
  probability plot of the fitted CDF), `nobs()`, `confint()` for all models.
* S3 methods are registered but no longer exported: use the generics (`summary(fit)`
  instead of `summary.Qest(fit)`).
* `Qest(..., Q = Qnorm)` no longer warns that `Qlm()` could be used.
* Weights are treated as sampling weights in the sandwich formula (the middle term
  is now the sum of w^2 g g', so that the standard errors do not depend on the scale
  of the weights).
* With left-truncated data, truncation is no longer ignored when all truncation
  times are smaller than the smallest observed time (in 1.0.2 this biased the
  estimates, e.g., with a common truncation point).
* `Qest()` with a user-defined `Q` warns, and uses a numerical Jacobian for the
  standard errors, if `Q` does not seem twice differentiable in `theta`.

## Dependencies

* 'pch', 'Hmisc' and 'segmented' code are no longer used; 'survival' and
  'matrixStats' moved from Depends to Imports.

## Documentation

* Documentation written with roxygen2, new vignette, tests (testthat), and a script
  (`data-raw/make_reference_1.0.2.R`) to compare with version 1.0.2.

# Qest 1.0.2

* Last version before the revision (CRAN, 2025-07-25).
