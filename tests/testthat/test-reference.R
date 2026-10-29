# Comparison with Qest 1.0.2 (tests/testthat/fixtures/ref_1.0.2.rds, created by
# data-raw/make_reference_1.0.2.R). Differences are expected where 1.0.2 was
# inaccurate (see NEWS.md): the tolerances are in units of standard errors.

refit <- function(cs, name) {
  r <- cs[[name]]
  env <- list2env(list(Qlin = Qlin, Qn = Qn, Qg = Qg, .refdata = r$data), parent = environment())
  call <- r$call
  if (identical(call[[1]], quote(suppressWarnings))) call <- call[[2]]
  call$data <- quote(.refdata)
  eval(call, env)
}

test_that("results agree with Qest 1.0.2", {
  skip_if_not(file.exists(ref_file()))
  cs <- ref_cases()
  tol <- c(norm_user = 0.03, norm_family = 1e-4, qlm = 1e-4, qlm_wtrunc = 1e-4,
           gamma_family = 0.03, gamma_user = 0.03, pois_family = 0.2, unif_min = 1e-3, unif_nomin = 1e-3,
           cens_family = 1e-3, cens_user = 0.05, trunc_user = 0.15, cox = 0.25)
  se_tol <- c(norm_family = 1e-4, qlm = 1e-4, qlm_wtrunc = 1e-4, gamma_family = 0.01, unif_min = 1e-3,
              unif_nomin = 1e-3)
  for (nm in intersect(names(tol), names(cs))) {
    fit <- suppressWarnings(refit(cs, nm))
    se <- cs[[nm]]$std.errs
    expect_lt(dse(coef(fit), cs[[nm]]$coefficients, se), tol[[nm]], label = paste("case", nm))
    if (nm %in% names(se_tol)) {
      se_new <- if (!is.null(fit$std.errs)) fit$std.errs else sqrt(diag(fit$var))
      expect_lt(max(abs(se_new / se - 1)), se_tol[[nm]], label = paste("SE, case", nm))
    }
  }
})
