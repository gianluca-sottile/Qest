# The generic engine: the Jacobian is the derivative of the estimating
# equations, and user-defined Q agree with the exact families.

ns <- asNamespace("Qest")
num_jac <- function(f, th, h = 1e-5) sapply(seq_along(th), function(j) {
  e <- replace(0 * th, j, h); (f(th + e) - f(th - e)) / (2 * h)
})

for (type in c("u", "c", "ct")) {
  test_that(paste("analytical Jacobian of the generic engine, type", type), {
    cs <- ref_cases()
    dd <- switch(type, u = cs$cens_family$data, c = cs$cens_family$data, ct = cs$trunc_user$data)
    n <- nrow(dd)
    ctx <- list(Q = Qn, data = dd, y = dd$y, z = if (type == "ct") dd$z else rep(-Inf, n),
                d = if (type == "u") rep(1, n) else dd$d, w = rep(1, n), type = type,
                tau = ns$buildTau(99, NULL, n, list()), h = 1e-4)
    th <- c(0.1, 1, 2)
    G <- ns$generic_ee(th, ctx, TRUE)
    Jn <- num_jac(function(t) ns$generic_ee(t, ctx)$g, th)
    expect_lt(max(abs(G$J - Jn)) / max(abs(Jn)), 1e-6)
  })
}

test_that("user-defined Q agrees with Qnorm (complete, censored, truncated data)", {
  cs <- ref_cases()
  d1 <- cs$norm_user$data
  m1 <- Qest(y ~ x1 + x2, Q = Qlin, start = rep(NA, 4), data = d1)
  m2 <- Qest(y ~ x1 + x2, Q = Qnorm, data = d1)
  expect_true(m1$converged)
  expect_lt(dse(coef(m1), coef(m2), m2$std.errs), 0.005)
  expect_equal(unname(m1$std.errs), unname(m2$std.errs), tolerance = 0.01)
  d5 <- cs$cens_family$data
  m1 <- Qest(Surv(y, d) ~ x, Q = Qn, start = rep(NA, 3), data = d5)
  m2 <- Qest(Surv(y, d) ~ x, Q = Qnorm, data = d5)
  expect_lt(dse(coef(m1), coef(m2), m2$std.errs), 0.005)
  expect_equal(unname(m1$std.errs), unname(m2$std.errs), tolerance = 0.01)
  d6 <- cs$trunc_user$data
  m1 <- Qest(Surv(z, y, d) ~ x, Q = Qn, start = rep(NA, 3), data = d6)
  m2 <- Qest(Surv(z, y, d) ~ x, Q = Qnorm, data = d6)
  expect_lt(dse(coef(m1), coef(m2), m2$std.errs), 0.005)
  expect_equal(unname(m1$std.errs), unname(m2$std.errs), tolerance = 0.01)
})

test_that("user-defined Q agrees with Qgamma", {
  d2 <- ref_cases()$gamma_user$data
  m2 <- Qest(y ~ x, Q = Qgamma, data = d2)
  m1 <- Qest(y ~ x, Q = Qg, start = coef(m2), data = d2)
  expect_lt(dse(coef(m1), coef(m2), m2$std.errs), 0.005)
  expect_equal(unname(m1$std.errs), unname(m2$std.errs), tolerance = 0.01)
})

test_that("the estimating equations of Qgamma are close to the exact ones", {
  d2 <- ref_cases()$gamma_user$data
  m <- Qest(y ~ x, Q = Qgamma, data = d2)
  th <- coef(m); a <- exp(th[1]); b <- exp(th[2] + th[3] * d2$x)
  # equation for the log-scale: int_0^1 A(t) dt - A(F(y)), with A(t) = b * a * pgamma(qgamma(t, a), a + 1)
  Aq <- function(t) a * pgamma(qgamma(t, a), a + 1)
  Aq1 <- integrate(Aq, 0, 1, rel.tol = 1e-10)$value
  g <- mean(b * (Aq1 - Aq(pgamma(d2$y / b, a))))
  expect_lt(abs(g), 1e-3 * mean(b))
})

test_that("Qnorm and Qlm give the same estimates", {
  d1 <- ref_cases()$norm_user$data
  m1 <- Qest(y ~ x1 + x2, Q = Qnorm, data = d1)
  m2 <- Qlm(y ~ x1 + x2, data = d1)
  expect_equal(unname(coef(m1)[-1]), unname(coef(m2)), tolerance = 1e-6)
  expect_equal(unname(m1$std.errs[-1]), unname(m2$std.errs), tolerance = 1e-5)
  expect_equal(unname(exp(coef(m1)[1])), as.numeric(m2$sigma), tolerance = 1e-6)
})

test_that("the starting points are estimated when 'start' contains NA", {
  set.seed(3)
  dd <- data.frame(x = runif(150))
  dd$y <- rgamma(150, shape = 3, scale = exp(0.5 + dd$x))
  m <- Qest(y ~ x, Q = Qg, start = c(NA, NA, NA), data = dd)
  expect_true(m$converged)
  expect_lt(max(abs(coef(m) - c(log(3), 0.5, 1)) / m$std.errs), 4)
})
