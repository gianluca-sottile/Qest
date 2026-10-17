test_that("the Frisch-Newton solver finds the quantile regression solution", {
  set.seed(1)
  n <- 300
  x <- cbind(1, rnorm(n), runif(n))
  y <- c(x %*% c(1, 2, 3)) + rt(n, 3)
  for (tau in c(0.1, 0.5, 0.9)) {
    f <- asNamespace("Qest")$rq_fnb(x, y, tau)
    obj <- function(b) sum((y - x %*% b) * (tau - (y - x %*% b < 0)))
    o <- optim(f$coefficients, obj, control = list(reltol = 1e-14, maxit = 5000))
    expect_lte(obj(f$coefficients), o$value + 1e-6)
    expect_lt(abs(mean(f$residuals < -1e-8) - tau), 0.02)
  }
})

test_that("wtrunc", {
  tau <- c(0.01, 0.5, 0.99)
  expect_equal(wtrunc(tau, 0.05, 0.05), c(0, 1, 0))
  expect_equal(wtrunc(tau, 0, 0.05), c(1, 1, 0))
  expect_equal(wtrunc(tau, 0.05, 0), c(0, 1, 1))
  expect_true(all(wtrunc(seq(0, 1, 0.01), smooth = TRUE) <= 1))
})

test_that("invQ inverts Q", {
  set.seed(1)
  dd <- data.frame(x = runif(20))
  u <- runif(20)
  y <- Qn(c(0, 1, 2), u, dd)
  expect_equal(invQ(Qn, c(0, 1, 2), y, dd, n.it = 30), u, tolerance = 1e-6)
})
