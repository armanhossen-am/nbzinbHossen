# Negative binomial regression by IRLS

Fits the NB2 regression model \\(\log \mu_i = x_i^\top \beta + o_i\\),
\\(\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta\\), by alternating (a)
iteratively reweighted least squares for \\(\beta\\) with \\(\theta\\)
fixed and (b) Newton-Raphson for \\(\log\theta\\) with \\(\beta\\)
fixed, until the log-likelihood changes by less than `control$epsilon`.
No general-purpose optimizer is used. The results agree with
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Usage

``` r
nbreg(
  formula,
  data,
  subset,
  na.action,
  weights,
  offset,
  control = nbreg_control(...),
  model = TRUE,
  ...
)
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`. Offsets can be given in the
  formula as `offset(log(exposure))`.

- data:

  A data frame containing the variables in the model.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values (default from `getOption("na.action")`).

- weights:

  Optional non-negative prior weights \\(a_i\\).

- offset:

  Optional offset vector (added to the linear predictor).

- control:

  A list from
  [`nbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg_control.md).

- model:

  Logical; keep the model frame in the result?

- ...:

  Arguments passed to
  [`nbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg_control.md)
  when `control` is not given.

## Value

An object of class `"nbreg"`, a list with components including
`coefficients`, `vcov` (the covariance of \\(\hat\beta\\), \\((X^\top
\hat W X)^{-1}\\), treating \\(\theta\\) as known), `theta`, `SE.theta`,
`loglik`, `iter`, `converged`, `fitted.values`, `linear.predictors`,
`y`, `x`, `weights`, `offset`, `df.residual` and `call`.

## See also

[`zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md),
[`summary.nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg-methods.md),
[`predict.nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg-methods.md)

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit)
#> 
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Deviance residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.1678 -1.3617 -0.2806  0.4476  3.4524 
#> 
#> Coefficients (NB2, log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .  
#> femWomen    -0.216418   0.072636  -2.979 0.002887 ** 
#> marMarried   0.150489   0.082097   1.833 0.066791 .  
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326    
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta: 2.264   Std. Err.: 0.2712 
#> Log-likelihood: -1560.96 on 7 df,  AIC: 3135.92 
#> Iterations (IRLS / Newton-Raphson alternations): 4 
#> 
fit$theta
#> [1] 2.264388
AIC(fit)
#> [1] 3135.917
```
