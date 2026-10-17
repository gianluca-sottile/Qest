#' Qest: Quantile-Based Estimators
#'
#' Quantile-based estimators (Q-estimators) fit a parametric model by
#' minimizing the integrated loss of quantile regression,
#' \deqn{\sum_i \int_0^1 w(\tau)\rho_\tau(y_i - Q(\tau \mid x_i, \theta))\,d\tau,}
#' where \eqn{Q(\tau \mid x, \theta)} is the quantile function of the model
#' (Sottile and Frumento, 2022). Q-estimators are consistent and asymptotically
#' normal, are usually more robust than maximum likelihood estimators, and allow
#' for censored and truncated data.
#'
#' The main functions are:
#' * [Qest()]: models defined by a user-written quantile function, or by one of
#'   the built-in families [Qnorm], [Qgamma], [Qpois] and [Qunif];
#' * [Qlm()]: linear models with normal errors (closed-form estimating equations);
#' * [Qcoxph()]: proportional hazards models with a piecewise linear baseline
#'   cumulative hazard.
#'
#' The weighting function \eqn{w(\tau)} (e.g., [wtrunc()]) can be used to
#' increase the robustness of the estimators.
#'
#' @references Sottile, G. and Frumento, P. (2022). Robust estimation and
#'   regression with parametric quantile functions. *Computational Statistics
#'   and Data Analysis*, 171, 107471. \doi{10.1016/j.csda.2022.107471}
#' @import stats
#' @importFrom survival survfit
#' @keywords internal
"_PACKAGE"
