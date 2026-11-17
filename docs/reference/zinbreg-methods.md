# Methods for zinbreg objects

Standard S3 methods for fitted zero-inflated negative binomial
regressions returned by
[`zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'zinbreg'
predict(
  object,
  newdata,
  type = c("response", "count", "zero", "prob"),
  at = NULL,
  ...
)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently unused).

- object, x:

  An object of class `"zinbreg"` (or `"summary.zinbreg"` for the print
  method of the summary).

- signif.stars:

  Logical; show significance stars?

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (count and
  zero coefficients, named `count_*` and `zero_*` as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)),
  `"count"` or `"zero"`.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"pearson"` (default) or `"response"`. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (marginal mean \\((1-\pi)\mu\\)), `"count"` (NB mean \\(\mu\\)),
  `"zero"` (structural-zero probability \\(\pi\\)) or `"prob"` (matrix
  of \\(P(Y = k)\\) for `k` in `at`).

- newdata:

  Optional data frame in which to look for covariates.

- at:

  For `type = "prob"`: the counts at which to evaluate the probabilities
  (default `0:max(y)`).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`"logLik"` with `df` = number of count and zero coefficients

- 1 (for \\(\theta\\));
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) numeric vectors (a
  matrix for `type = "prob"`);
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
  observations with positive weight;
  [`summary()`](https://rdrr.io/r/base/summary.html) an object of class
  `"summary.zinbreg"` with coefficient tables for both model parts.

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
coef(fit, model = "zero")
#> (Intercept)        ment 
#>  -0.7949062  -0.6193810 
AIC(fit)
#> [1] 3123.505
nd <- data.frame(fem = "Men", kid5 = 0, ment = c(0, 10, 40))
predict(fit, newdata = nd, type = "response")
#>        1        2        3 
#> 1.128592 2.088116 4.339573 
predict(fit, newdata = nd, type = "zero")
#>            1            2            3 
#> 3.111162e-01 9.213804e-04 7.852850e-12 
# expected number of zeros
sum(predict(fit, type = "prob", at = 0))
#> [1] 286.3748
```
