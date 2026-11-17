# Methods for nbreg objects

Standard S3 methods for fitted negative binomial regressions returned by
[`nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("deviance", "pearson", "response"), ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("link", "response"), ...)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently unused).

- object, x:

  An object of class `"nbreg"` (or `"summary.nbreg"` for the print
  method of the summary).

- signif.stars:

  Logical; show significance stars?

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"deviance"` (default), `"pearson"` or `"response"`. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"link"` (linear
  predictor \\(\eta\\)) or `"response"` (mean \\(\mu\\)).

- newdata:

  Optional data frame in which to look for covariates (and offset
  variables) for prediction. If omitted, the fitted values are used.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`"logLik"` with `df` = number of coefficients + 1 (for \\(\theta\\)), so
that [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) numeric vectors;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of observations
with positive weight; [`summary()`](https://rdrr.io/r/base/summary.html)
an object of class `"summary.nbreg"` with a coefficient table of
estimates, standard errors, z values and p-values.

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- nbreg(art ~ fem + kid5 + ment, data = bioChemists)
fit
#> 
#> Call:
#> nbreg(formula = art ~ fem + kid5 + ment, data = bioChemists)
#> 
#> Negative binomial (NB2) regression, log link, fitted by IRLS
#> 
#> Coefficients:
#> (Intercept)     femWomen         kid5         ment  
#>     0.39102     -0.23270     -0.13775      0.02937  
#> 
#> Theta: 2.241   (SE 0.2668 )
#> Log-likelihood: -1562.66 on 5 df
#> Iterations: 4 
#> 
coef(fit)
#> (Intercept)    femWomen        kid5        ment 
#>  0.39102467 -0.23269704 -0.13775445  0.02936915 
sqrt(diag(vcov(fit)))
#> (Intercept)    femWomen        kid5        ment 
#> 0.064526913 0.072183539 0.048153108 0.003116131 
logLik(fit)
#> 'log Lik.' -1562.664 (df=5)
BIC(fit)
#> [1] 3159.423
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -1.0015086 -0.9276857 -0.9276857 -0.9296379 -1.0884583 -0.8147337 
predict(fit, newdata = data.frame(fem = "Women", kid5 = 0, ment = 10),
        type = "response")
#>        1 
#> 1.571482 
```
