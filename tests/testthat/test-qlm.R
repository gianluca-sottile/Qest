test_that("Qlm: estimates, methods and weights", {
  d1 <- ref_cases()$norm_user$data
  m <- Qlm(y ~ x1 + x2, data = d1)
  expect_s3_class(m, "Qlm")
  expect_true(m$converged)
  expect_equal(unname(m$residuals + m$fitted.values), d1$y)
  expect_equal(predict(m), m$fitted.values)
  nd <- data.frame(x1 = c(0, 1), x2 = c(1, 2))
  p <- predict(m, nd, se.fit = TRUE)
  expect_equal(unname(p$fit), unname(c(cbind(1, as.matrix(nd)) %*% coef(m))))
  expect_equal(unname(predict(m, nd, type = "quantile", p = 0.5)[, 1]), unname(p$fit))
  expect_equal(vcov(m), m$covar)
  expect_equal(nobs(m), 100L)
  expect_output(print(summary(m)), "sigma")
  ci <- confint(m)
  expect_true(all(ci[, 1] < coef(m) & ci[, 2] > coef(m)))
  # integer weights are equivalent to replicated observations
  w <- rep(1:2, 50)
  mw <- Qlm(y ~ x1 + x2, data = d1, weights = w)
  mr <- Qlm(y ~ x1 + x2, data = d1[rep(1:100, w), ])
  expect_equal(coef(mw), coef(mr), tolerance = 1e-6)
})

test_that("Qlm is robust with truncated weights", {
  set.seed(5)
  x <- runif(200)
  y <- 1 + 2 * x + rnorm(200)
  y[1:10] <- y[1:10] + 30
  m0 <- Qlm(y ~ x)
  m1 <- Qlm(y ~ x, wtau = wtrunc, delta.left = 0.1, delta.right = 0.1)
  expect_lt(abs(coef(m1)[1] - 1), abs(coef(lm(y ~ x))[1] - 1))
  expect_true(m1$converged)
})

test_that("Qlm: start, errors, collinearity", {
  d1 <- ref_cases()$norm_user$data
  m <- Qlm(y ~ x1 + x2, data = d1)
  m2 <- Qlm(y ~ x1 + x2, data = d1, start = coef(m) + 1)
  expect_equal(coef(m), coef(m2), tolerance = 1e-6)
  expect_error(Qlm(y ~ x1, data = d1, start = 1:3), "length")
  expect_error(Qlm(y ~ 0, data = d1), "no covariates")
  expect_error(Qlm(y ~ x1, data = d1, foo = 1), "only allowed")
  d1$x3 <- d1$x1
  expect_warning(m3 <- Qlm(y ~ x1 + x3, data = d1), "collinearity")
  expect_true(is.na(coef(m3)["x3"]))
})
