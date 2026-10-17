
<!-- README.md is generated from README.Rmd. Please edit that file -->


# Qest

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/Qest)](https://CRAN.R-project.org/package=Qest)
[![R-CMD-check](https://github.com/gianluca-sottile/Qest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gianluca-sottile/Qest/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

'Qest' implements the quantile-based estimators (Q-estimators) of Sottile and
Frumento (2022): any parametric model is fitted through its quantile function by
minimizing the integrated loss of quantile regression. Q-estimators are usually
more robust than maximum likelihood, and allow for censored and truncated data.

* `Qest()`: user-defined quantile functions, or the families `Qnorm`, `Qgamma`,
  `Qpois` (continuous Poisson) and `Qunif`;
* `Qlm()`: linear models with normal errors;
* `Qcoxph()`: proportional hazards models with a piecewise linear baseline hazard.

## Installation

```r
install.packages("Qest")
# development version
remotes::install_github("gianluca-sottile/Qest")
```

## Example

``` r
library(Qest)
set.seed(1)
n <- 200
dat <- data.frame(x = runif(n))
dat$y <- rgamma(n, shape = 2, scale = exp(1 + dat$x))
Qg <- function(theta, tau, data) {
  qgamma(tau, shape = exp(theta[1]), scale = exp(theta[2] + theta[3] * data$x))
}
m1 <- Qest(y ~ x, Q = Qg, start = c(NA, NA, NA), data = dat)
m2 <- Qest(y ~ x, Q = Qgamma, data = dat)
cbind(user = coef(m1), family = coef(m2))
#>             user    family
#> Theta1 0.6636152 0.6635124
#> Theta2 0.9955430 0.9956836
#> Theta3 0.9858171 0.9858114
summary(m2)
#> 
#> Call: Qest(formula = y ~ x, Q = Qgamma, data = dat)
#> 
#> n. of observations: 200 
#> n. of parameters: 3 
#> n. of iterations: 4 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> log(shape)    0.6635     0.1084   6.120 9.37e-10 ***
#> (Intercept)   0.9957     0.1650   6.034 1.59e-09 ***
#> x             0.9858     0.1933   5.099 3.41e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Minimized loss function: 331.7161
#> 
```

## Reference

Sottile, G. and Frumento, P. (2022). Robust estimation and regression with
parametric quantile functions. *Computational Statistics and Data Analysis*, 171,
107471. <https://doi.org/10.1016/j.csda.2022.107471>
