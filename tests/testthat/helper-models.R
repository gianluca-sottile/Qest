# Quantile functions and data used in several tests.
library(survival)

Qlin <- function(theta, tau, data) {
  X <- model.matrix(~ x1 + x2, data = data)
  qnorm(tau, c(X %*% theta[-1]), exp(theta[1]))
}
Qn <- function(theta, tau, data) qnorm(tau, theta[2] + theta[3] * data$x, exp(theta[1]))
Qg <- function(theta, tau, data) qgamma(tau, shape = exp(theta[1]), scale = exp(theta[2] + theta[3] * data$x))

ref_file <- function() test_path("fixtures", "ref_1.0.2.rds")
ref_cases <- function() readRDS(ref_file())$cases

# relative difference in units of standard errors
dse <- function(a, b, se) max(abs((a - b) / se))

sim_cox <- function(n = 200, seed = 7) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 3), x2 = rbinom(n, 1, 0.5))
  t <- rweibull(n, shape = 2, scale = (1 / exp(2 + 2 * d$x - d$x2))^(1 / 2))
  cens <- runif(n, 0, 1)
  d$y <- pmin(t, cens); d$d <- as.numeric(t <= cens)
  d
}
