# Grid of quantile orders used for numerical integration.
#
# tau[1] = omicron, tau[ntau] = 1 - omicron, and ntau - 2 points in between that
# are equally spaced on the scale of a Beta(2,2) quantile function (denser in the
# tails). dtau are half steps (trapezoids), as in Qest <= 1.0.2, so that the code
# of the families can be used unchanged. tau[2] = 1.1*omicron (see findAp).
buildTau <- function(ntau, wtau = NULL, nobs, wtauoptions = NULL, omicron = 5e-4) {
  if (!is.numeric(ntau) || length(ntau) != 1 || ntau < 9) stop("'ntau' must be a number >= 9", call. = FALSE)
  ntau <- round(ntau)
  tauL <- omicron; tauR <- 1 - omicron
  tau <- stats::pbeta(seq.int(stats::qbeta(omicron * 1.1, 2, 2), stats::qbeta(1 - omicron * 1.1, 2, 2),
                              length.out = ntau - 2), 2, 2)
  tau <- c(tauL, tau, tauR)
  dtau <- 0.5 * (c(tau[1], tau[2:ntau] - tau[1:(ntau - 1)]))

  wfun <- wtau
  wfit <- !is.null(wfun)
  if (!wfit) wfun <- function(tau, ...) rep(1, length(tau))
  wt <- callwtau(wfun, tau, wtauoptions)
  w01 <- callwtau(wfun, 0:1, wtauoptions)
  if (wfit) {
    if (length(wt) != ntau || anyNA(wt) || any(wt < 0) || any(!is.finite(c(wt, w01))))
      stop("'wtau(tau)' must return a non-negative, finite value for 0 <= tau <= 1", call. = FALSE)
    if (all(wt == 0)) stop("'wtau(tau)' is zero everywhere", call. = FALSE)
    wt <- pmax(wt, 1e-12)
  }

  list(tau = tau, dtau = dtau, ntau = ntau, TAU = t(matrix(tau, ntau, nobs)),
       tau1 = c(0, tau, 1), dtau1 = c(dtau, 1 - tau[ntau]), ntau1 = ntau + 2,
       dtau_long = rep(dtau, rep.int(nobs, ntau)), wtau_long = rep(wt, rep.int(nobs, ntau)),
       wtau = wt, wfun = wfun, wfit = wfit, tauL = tauL, tauR = tauR, opt = wtauoptions,
       omicron = omicron)
}
