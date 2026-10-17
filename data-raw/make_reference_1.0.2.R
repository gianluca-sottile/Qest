# Reference results of Qest 1.0.2 (CRAN), used by tests/testthat/test-reference.R.
#
# Run with Qest 1.0.2 installed (NOT the development version):
#   install.packages("Qest")     # or remotes::install_version("Qest", "1.0.2")
#   source("data-raw/make_reference_1.0.2.R")
# The fixture is saved in tests/testthat/fixtures/ref_1.0.2.rds.
# Version 2.0.0 changes the numerical derivatives of the generic engine and the
# starting points, so the comparison is made with a tolerance (see the test file).

library(Qest)
stopifnot(packageVersion("Qest") == "1.0.2")

cases <- list()
run <- function(name, expr, data) {
  fit <- try(eval(expr), silent = TRUE)
  if (inherits(fit, "try-error")) { message(name, ": ", fit); return(invisible()) }
  cases[[name]] <<- list(call = expr, data = data, coefficients = coef(fit),
                        std.errs = if (!is.null(fit$std.errs)) fit$std.errs else sqrt(diag(fit$var)),
                        obj.function = as.numeric(fit$obj.function))
  message(name, " done")
}

# 1. Normal linear model: user-defined Q, Qnorm, Qlm
set.seed(1)
n <- 100
d1 <- data.frame(x1 = rnorm(n), x2 = runif(n, 0, 3))
d1$y <- 4 + d1$x1 + 2 * d1$x2 + exp(1) * qnorm(runif(n))
Qlin <- function(theta, tau, data) {
  X <- model.matrix(~ x1 + x2, data = data)
  qnorm(tau, X %*% theta[-1], exp(theta[1]))
}
run("norm_user", quote(Qest(y ~ x1 + x2, Q = Qlin, start = rep(NA, 4), data = d1)), d1)
run("norm_family", quote(suppressWarnings(Qest(y ~ x1 + x2, Q = Qnorm, data = d1))), d1)
run("qlm", quote(Qlm(y ~ x1 + x2, data = d1)), d1)
run("qlm_wtrunc", quote(Qlm(y ~ x1 + x2, data = d1, wtau = wtrunc, delta.left = 0.05, delta.right = 0.05)), d1)

# 2. Gamma
set.seed(2)
n <- 250
d2 <- data.frame(x = runif(n))
d2$y <- rgamma(n, shape = exp(1), scale = exp(1 + 2 * d2$x))
Qg <- function(theta, tau, data) {
  qgamma(tau, shape = exp(theta[1]), scale = exp(theta[2] + theta[3] * data$x))
}
run("gamma_family", quote(Qest(y ~ x, Q = Qgamma, data = d2)), d2)
run("gamma_user", quote(suppressWarnings(Qest(y ~ x, Q = Qg, start = rep(NA, 3), data = d2))), d2)

# 3. Poisson (jittered) and Uniform
set.seed(3)
d3 <- data.frame(x = runif(n))
d3$y <- rpois(n, exp(1 + 2 * d3$x)) + runif(n)
run("pois_family", quote(Qest(y ~ x, Q = Qpois, data = d3)), d3)
set.seed(4)
d4 <- data.frame(x = runif(n))
d4$y <- runif(n, 0, 1 + 2 * d4$x)
run("unif_min", quote(Qest(y ~ x, Q = Qunif(min = TRUE), data = d4)), d4)
run("unif_nomin", quote(Qest(y ~ x, Q = Qunif(min = FALSE), data = d4)), d4)

# 4. Censored and truncated normal
set.seed(5)
d5 <- data.frame(x = runif(n))
t5 <- rnorm(n, 1 + 2 * d5$x, 1); c5 <- rnorm(n, 3, 1)
d5$y <- pmin(t5, c5); d5$d <- as.numeric(t5 <= c5)
Qn <- function(theta, tau, data) qnorm(tau, theta[2] + theta[3] * data$x, exp(theta[1]))
run("cens_family", quote(Qest(Surv(y, d) ~ x, Q = Qnorm, data = d5)), d5)
run("cens_user", quote(Qest(Surv(y, d) ~ x, Q = Qn, start = rep(NA, 3), data = d5)), d5)
set.seed(6)
m <- 600
x6 <- runif(m); t6 <- rnorm(m, 1 + 2 * x6, 1); z6 <- rnorm(m, 0.5, 1); c6 <- rnorm(m, 4, 1)
keep <- t6 > z6 & c6 > z6
d6 <- data.frame(x = x6, z = z6, y = pmin(t6, c6), d = as.numeric(t6 <= c6))[keep, ][1:250, ]
run("trunc_user", quote(Qest(Surv(z, y, d) ~ x, Q = Qn, start = rep(NA, 3), data = d6)), d6)

# 5. Cox model (needs the real 'survival' package)
set.seed(7)
d7 <- data.frame(x = runif(200, 0, 3))
t7 <- rweibull(200, shape = 2, scale = (1 / exp(2 + 2 * d7$x))^(1 / 2))
c7 <- runif(200, 0, 1)
d7$y <- pmin(t7, c7); d7$d <- as.numeric(t7 <= c7)
run("cox", quote(Qcoxph(Surv(y, d) ~ x, data = d7)), d7)

ref <- list(version = as.character(packageVersion("Qest")), R = R.version.string,
            date = Sys.Date(), cases = cases)
saveRDS(ref, "tests/testthat/fixtures/ref_1.0.2.rds", version = 2)
