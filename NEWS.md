# nbzinbHossen 0.1.0

First release, submitted for the STA 651 final project.

## Model fitting

* `nbreg()` fits negative binomial (NB2) regression by IRLS for the
  coefficients alternated with Newton-Raphson on `log(theta)`. It agrees with
  `MASS::glm.nb()` to about 1e-7.
* `zinbreg()` fits zero-inflated negative binomial regression by an EM
  algorithm whose M-steps are IRLS fits, with two-part formulas
  `y ~ count | zero`. It agrees with `pscl::zeroinfl(dist = "negbin")` to
  about 1e-5.
* ZINB standard errors come from an analytic Hessian of the observed-data
  log-likelihood; `zinbreg_control(hessian = "numerical")` uses
  `numDeriv::hessian()` instead.
* No general-purpose optimizer (`optim()`, `nlm()`, ...) and no `glm()` /
  `lm()` fitting code is used.

## Methods

* `print()`, `summary()`, `coef()`, `vcov()`, `logLik()`, `fitted()`,
  `residuals()`, `nobs()` and `predict()` for both classes; `AIC()` and
  `BIC()` work through `logLik()`.
* `predict.zinbreg()` supports `type = "response"`, `"count"`, `"zero"` and
  `"prob"`, with `newdata`.

## Checks

* 82 tests compare both fitting functions with the reference implementations
  and cover offsets, weights, warnings and input errors.
* Benchmark and simulation scripts are in `inst/benchmark/`.
