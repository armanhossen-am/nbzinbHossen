# Zero-inflated negative binomial regression by EM

Fits the ZINB model \$\$P(Y_i = 0) = \pi_i + (1 - \pi_i) f\_{NB}(0;
\mu_i, \theta), \quad P(Y_i = y) = (1 - \pi_i) f\_{NB}(y; \mu_i,
\theta),\\ y \ge 1,\$\$ with \\(\log\mu_i = x_i^\top\beta + o_i\\) and
\\(\mathrm{logit}\\,\pi_i = z_i^\top\gamma\\), by an EM algorithm. The
E-step computes the posterior probability \\(\tau_i\\) that observation
\\(i\\) is a structural zero; the M-step fits a logistic regression with
fractional response \\(\tau_i\\) by IRLS (zero part) and a weighted NB
regression with prior weights \\(1 - \tau_i\\) by
[`nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md)'s
IRLS / Newton-Raphson algorithm (count part). Standard errors come from
the observed information of the observed-data log-likelihood in
\\((\beta, \gamma, \log\theta)\\). The results agree with
`pscl::zeroinfl(..., dist = "negbin")`.

## Usage

``` r
zinbreg(
  formula,
  data,
  subset,
  na.action,
  weights,
  offset,
  control = zinbreg_control(...),
  model = TRUE,
  ...
)
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before `|`, zero-model covariates after it. Without `|` the same
  covariates are used in both parts. Offsets (in the formula or via
  `offset`) apply to the count part only.

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
  [`zinbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg_control.md).

- model:

  Logical; keep the model frame in the result?

- ...:

  Arguments passed to
  [`zinbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg_control.md)
  when `control` is not given.

## Value

An object of class `"zinbreg"`, a list with components including
`coefficients` (a list with `count` and `zero`), `vcov` (covariance of
all count and zero coefficients), `vcov_full` (also including
`log(theta)`), `theta`, `SE.logtheta`, `SE.theta`, `loglik`, `iter`,
`converged`, `loglik_trace`, `fitted.values` (\\((1-\pi_i)\mu_i\\)),
`mu`, `pi`, `y`, `x`, `z`, `weights`, `offset`, `df.residual` and
`call`.

## See also

[`nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md),
[`summary.zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg-methods.md),
[`predict.zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg-methods.md)

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
               data = bioChemists)
summary(fit)
#> 
#> Call:
#> zinbreg(formula = art ~ fem + mar + kid5 + phd + ment | ment, 
#>     data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3041 -0.7685 -0.2632  0.4671  6.3765 
#> 
#> Count model coefficients (negbin with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.404005   0.141715   2.851  0.00436 ** 
#> femWomen    -0.211906   0.071922  -2.946  0.00322 ** 
#> marMarried   0.139462   0.081193   1.718  0.08586 .  
#> kid5        -0.167624   0.052457  -3.195  0.00140 ** 
#> phd          0.001963   0.035586   0.055  0.95601    
#> ment         0.024393   0.003518   6.934 4.08e-12 ***
#> Log(theta)   1.002892   0.142829   7.022 2.19e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.8067     0.3532  -2.284   0.0224 *
#> ment         -0.6095     0.2458  -2.480   0.0131 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta = 2.726   (SE 0.3894 )
#> Number of EM iterations: 139 
#> Log-likelihood: -1553.27 on 9 df,  AIC: 3124.54 
#> 
head(predict(fit, type = "zero"))
#>            1            2            3            4            5            6 
#> 6.223342e-03 1.138847e-02 1.138847e-02 6.690807e-02 5.852534e-08 1.165331e-01 
```
