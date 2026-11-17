# Changelog

## nbzinbHossen 0.1.0

First release, submitted for the STA 651 final project.

### Model fitting

- [`nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md)
  fits negative binomial (NB2) regression by IRLS for the coefficients
  alternated with Newton-Raphson on `log(theta)`. It agrees with
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) to about
  1e-7.
- [`zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md)
  fits zero-inflated negative binomial regression by an EM algorithm
  whose M-steps are IRLS fits, with two-part formulas
  `y ~ count | zero`. It agrees with `pscl::zeroinfl(dist = "negbin")`
  to about 1e-5.
- ZINB standard errors come from an analytic Hessian of the
  observed-data log-likelihood; `zinbreg_control(hessian = "numerical")`
  uses
  [`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
  instead.
- No general-purpose optimizer
  ([`optim()`](https://rdrr.io/r/stats/optim.html),
  [`nlm()`](https://rdrr.io/r/stats/nlm.html), …) and no
  [`glm()`](https://rdrr.io/r/stats/glm.html) /
  [`lm()`](https://rdrr.io/r/stats/lm.html) fitting code is used.

### Methods

- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) for both classes;
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) work through
  [`logLik()`](https://rdrr.io/r/stats/logLik.html).
- [`predict.zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg-methods.md)
  supports `type = "response"`, `"count"`, `"zero"` and `"prob"`, with
  `newdata`.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for `zinbreg`
  fits shows the observed-data log-likelihood at every EM iteration, or
  its distance to the final value on a log scale.

### Checks

- 87 tests compare both fitting functions with the reference
  implementations and cover offsets, weights, warnings and input errors.
- Benchmark and simulation scripts are in `inst/benchmark/`.
