# nbzinbHossen

<!-- badges: start -->
<!-- badges: end -->

**nbzinbHossen** fits

- **negative binomial (NB2) regression** with `nbreg()` (an alternative to
  `MASS::glm.nb()`), using iteratively reweighted least squares (IRLS) for
  the coefficients alternated with Newton–Raphson on log(θ); and
- **zero-inflated negative binomial (ZINB) regression** with `zinbreg()` (an
  alternative to `pscl::zeroinfl(..., dist = "negbin")`), using an EM
  algorithm whose M-steps are IRLS fits (a fractional-response logistic
  regression for the zero part and a weighted NB regression for the count
  part).

No general-purpose optimizer (`optim()`, `nlm()`, …) and no `glm()`/`lm()`
fitting code is used. ZINB standard errors come from an **analytic** Hessian of
the observed-data log-likelihood.

Website (function reference + project report):
<https://armanhossen-am.github.io/nbzinbHossen/>

## Installation

```r
# install.packages("devtools")
devtools::install_github("armanhossen-am/nbzinbHossen")
```

## Usage

```r
library(nbzinbHossen)
data("bioChemists", package = "pscl")

# Negative binomial regression
nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(nb)
nb$theta                     # 2.264, same as MASS::glm.nb()

# Zero-inflated NB: count covariates | zero-inflation covariates
zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                data = bioChemists)
summary(zinb)
AIC(nb, zinb)

# Predictions: marginal mean, NB mean, and structural-zero probability
nd <- data.frame(fem = "Women", mar = "Married", kid5 = 1, phd = 3, ment = 10)
predict(zinb, newdata = nd, type = "response")
predict(zinb, newdata = nd, type = "count")
predict(zinb, newdata = nd, type = "zero")
```

Both fitting functions return S3 objects with `print`, `summary`, `coef`,
`vcov`, `logLik` (so `AIC()`/`BIC()` work), `fitted`, `residuals`, `nobs` and
`predict` methods. Use `nbreg_control()` / `zinbreg_control()` to change
tolerances, iteration limits, or to print the iteration history
(`trace = TRUE`).

## Agreement with the reference implementations

On `pscl::bioChemists` and `AER::NMES1988` the maximum absolute differences
are about 1e-7 from `MASS::glm.nb()` for NB, and below 3e-5 from
`pscl::zeroinfl()` for ZINB. The tolerances required by the project are
1e-5 and 1e-3. The full tables are in `inst/benchmark/` and in the report.

## R CMD check

`devtools::check()` gives 0 errors, 0 warnings and 0 notes.
