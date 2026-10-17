test_that("the families recover the parameters", {
  set.seed(11)
  n <- 1000
  x <- runif(n)
  eta <- 1 + 2 * x
  dd <- data.frame(x = x)
  dd$yn <- rnorm(n, eta, exp(0.5))
  dd$yg <- rgamma(n, shape = exp(1), scale = exp(eta))
  dd$yp <- rpois(n, exp(eta)) + runif(n)
  dd$yu <- runif(n, x, eta)
  m <- Qest(yn ~ x, Qnorm, data = dd)
  expect_lt(max(abs(coef(m) - c(0.5, 1, 2)) / m$std.errs), 4)
  m <- Qest(yg ~ x, Qgamma, data = dd)
  expect_lt(max(abs(coef(m) - c(1, 1, 2)) / m$std.errs), 4)
  m <- Qest(yp ~ x, Qpois, data = dd)
  expect_true(m$converged)
  expect_lt(max(abs(coef(m) - c(1, 2)) / m$std.errs), 4)
  m <- Qest(yu ~ x, Qunif(min = TRUE), data = dd)
  expect_lt(max(abs(coef(m) - c(0, 1, 1, 2)) / m$std.errs), 4)
  expect_named(coef(m), c("min: (Intercept)", "min: x", "max: (Intercept)", "max: x"))
})

test_that("Qfamily can be given as a name, a function or a call", {
  dd <- ref_cases()$gamma_user$data
  m1 <- Qest(y ~ x, "Qgamma", data = dd)
  m2 <- Qest(y ~ x, Qgamma, data = dd)
  m3 <- Qest(y ~ x, Qgamma(), data = dd)
  expect_equal(coef(m1), coef(m2))
  expect_equal(coef(m1), coef(m3))
  expect_output(print(Qunif(min = FALSE)), "Qunif")
})

test_that("the continuous Poisson family agrees with the exact quantile function", {
  skip_on_cran()
  dd <- ref_cases()$pois_family$data
  Qpc <- function(theta, tau, data) {
    lam <- exp(theta[1] + theta[2] * data$x)
    f <- function(tt, lam) {
      lo <- qpois(tt, lam); hi <- lo + 1
      for (i in 1:40) { mid <- (lo + hi) / 2; up <- pgamma(lam, mid, lower.tail = FALSE) < tt; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      (lo + hi) / 2
    }
    if (is.matrix(tau)) matrix(f(c(tau), rep(lam, ncol(tau))), nrow(tau)) else f(tau, lam)
  }
  m1 <- Qest(y ~ x, Q = Qpois, data = dd)
  m2 <- Qest(y ~ x, Q = Qpc, start = coef(m1), data = dd)
  expect_lt(dse(coef(m1), coef(m2), m1$std.errs), 0.01)
  expect_equal(unname(m1$std.errs), unname(m2$std.errs), tolerance = 0.01)
})

test_that("Poisson offsets", {
  set.seed(2)
  dd <- data.frame(x = runif(300), t = runif(300, 1, 3))
  dd$y <- rpois(300, dd$t * exp(1 + dd$x)) + runif(300)
  m1 <- Qest(y ~ x + offset(log(t)), Qpois, data = dd)
  m2 <- Qest(y ~ x, Qpois(offset = log(dd$t)), data = dd)
  expect_equal(coef(m1), coef(m2))
  expect_error(Qest(y ~ x + offset(log(t)), Qgamma, data = dd), "offsets")
})

test_that("unsupported data and models give errors", {
  dd <- ref_cases()$cens_family$data
  expect_error(Qest(Surv(y, d) ~ x, Qpois, data = dd), "not supported")
  expect_error(Qest(Surv(y, d) ~ x, Qunif, data = dd), "not supported")
  expect_error(Qest(y ~ x, "Qfoo", data = dd), "not recognized")
})

test_that("collinear covariates are dropped with a warning", {
  dd <- ref_cases()$gamma_user$data
  dd$x2 <- 2 * dd$x
  expect_warning(m <- Qest(y ~ x + x2, Qgamma, data = dd), "collinearity")
  expect_true(is.na(coef(m)["x2"]))
})

test_that("truncation is not ignored when all truncation times are below the observed times", {
  set.seed(9)
  y <- rnorm(5000, 1, 1); y <- y[y > 0][1:1000]
  m <- Qest(Surv(rep(0, 1000), y, rep(1, 1000)) ~ 1, Qnorm)
  expect_equal(m$internal$type, "ct")
  expect_lt(max(abs(coef(m) - c(0, 1)) / m$std.errs), 4)
})

test_that("predictions with Poisson offsets", {
  set.seed(2)
  dd <- data.frame(x = runif(300), t = runif(300, 1, 3))
  dd$y <- rpois(300, dd$t * exp(1 + dd$x)) + runif(300)
  m1 <- Qest(y ~ x + offset(log(t)), Qpois, data = dd)
  lam <- dd$t[1:3] * exp(coef(m1)[1] + coef(m1)[2] * dd$x[1:3])
  p1 <- predict(m1, newdata = dd[1:3, ], p = 0.5)
  p0 <- predict(m1, p = 0.5)[1:3, ]
  expect_equal(unname(p1[, 1]), unname(p0))
  expect_equal(unname(p0), asNamespace("Qest")$qpois_cont(0.5, lam), tolerance = 1e-6)
  m2 <- Qest(y ~ x, Qpois(offset = log(dd$t)), data = dd)
  expect_error(predict(m2, newdata = dd[1:3, ]), "offset")
})

test_that("standard errors do not depend on the scale of the weights", {
  dd <- ref_cases()$gamma_user$data
  w <- runif(nrow(dd), 0.5, 2)
  m1 <- Qest(y ~ x, Qgamma, data = dd, weights = w)
  m2 <- Qest(y ~ x, Qgamma, data = dd, weights = 10 * w)
  expect_equal(coef(m1), coef(m2), tolerance = 1e-8)
  expect_equal(m1$std.errs, m2$std.errs, tolerance = 1e-6)
})
