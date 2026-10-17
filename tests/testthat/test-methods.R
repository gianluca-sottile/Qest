test_that("Qest methods", {
  d1 <- ref_cases()$norm_user$data
  m1 <- Qest(y ~ x1 + x2, Q = Qlin, start = rep(NA, 4), data = d1)
  m2 <- Qest(y ~ x1 + x2, Q = Qnorm, data = d1)
  expect_output(print(m1), "Coefficients")
  expect_output(print(summary(m1, covar = TRUE)), "Covariance")
  expect_equal(vcov(m1), m1$covar)
  expect_equal(nobs(m1), 100L)
  nd <- d1[1:3, ]
  p1 <- predict(m1, nd, p = c(0.1, 0.5)); p2 <- predict(m2, nd, p = c(0.1, 0.5))
  expect_equal(dim(p1), c(3L, 2L))
  expect_equal(p1, p2, tolerance = 1e-3)
  expect_equal(unname(predict(m2, nd, type = "CDF")), unname(m2$CDF[1:3]), tolerance = 1e-6)
  expect_equal(unname(predict(m1, nd, type = "CDF")), unname(m1$CDF[1:3]), tolerance = 1e-3)
  pdf(NULL); on.exit(grDevices::dev.off())
  km <- plot(m1)
  expect_true(all(diff(km$F) >= 0))
})

test_that("rows of 'data' with missing values or zero weights are not passed to Q", {
  d1 <- ref_cases()$norm_user$data
  d1$x1[3] <- NA
  w <- rep(1, 100); w[5] <- 0
  expect_warning(m <- Qest(y ~ x1 + x2, Q = Qlin, start = rep(NA, 4), data = d1, weights = w), "zero weight")
  expect_equal(nrow(m$data), 98L)
  expect_length(m$CDF, 98)
})

test_that("informative errors", {
  d1 <- ref_cases()$norm_user$data
  expect_error(Qest(y ~ x1, Q = Qlin, data = d1), "start")
  expect_error(Qest(y ~ x1, Q = Qlin, start = rep(NA, 4), data = d1, foo = 1), "unused")
  expect_error(Qest(y ~ x1, Q = function(a, b, c) a, start = 1, data = d1), "arguments")
  expect_error(Qest.control(safeit = 10), "removed")
  expect_error(Qest.control(alpha0 = 1), "removed")
  expect_error(Qest(y ~ x1, Q = Qnorm, data = d1, weights = rep(-1, 100)), "non-negative")
})

test_that("wtau arguments are passed through ...", {
  set.seed(4)
  y <- rgamma(200, shape = 3, scale = 2)
  m1 <- Qest(y ~ 1, Q = Qgamma, wtau = function(tau, h) dnorm((tau - 0.5) / h), h = 0.2)
  m2 <- Qest(y ~ 1, Q = Qgamma, wtau = function(tau, h) dnorm((tau - 0.5) / h), h = 0.4)
  expect_false(isTRUE(all.equal(coef(m1), coef(m2))))
  expect_error(Qest(y ~ 1, Q = Qgamma, wtau = function(x) x), "first argument")
})
