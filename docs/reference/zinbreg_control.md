# Control parameters for zinbreg()

Control parameters for zinbreg()

## Usage

``` r
zinbreg_control(
  epsilon = 1e-10,
  maxit = 5000,
  nb_control = nbreg_control(epsilon = 1e-10),
  maxit_logit = 100,
  epsilon_logit = 1e-12,
  hessian = c("analytic", "numerical"),
  trace = FALSE
)
```

## Arguments

- epsilon:

  Convergence tolerance for the EM algorithm: iteration stops when the
  observed-data log-likelihood changes by less than `epsilon`.

- maxit:

  Maximum number of EM iterations.

- nb_control:

  Control list from
  [`nbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg_control.md)
  used for the initial NB fit and for the count-part M-step.

- maxit_logit, epsilon_logit:

  Maximum iterations and relative tolerance of the IRLS logistic
  regression in the zero-part M-step.

- hessian:

  How to compute the observed information for the standard errors:
  `"analytic"` (default, closed-form second derivatives) or
  `"numerical"`
  ([`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
  of the observed-data log-likelihood; requires the numDeriv package).

- trace:

  Logical; print the observed-data log-likelihood at each EM iteration?

## Value

A list of control parameters.

## Examples

``` r
zinbreg_control(epsilon = 1e-10, trace = TRUE)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
#> [1] 5000
#> 
#> $nb_control
#> $nb_control$epsilon
#> [1] 1e-10
#> 
#> $nb_control$maxit
#> [1] 100
#> 
#> $nb_control$epsilon_inner
#> [1] 1e-10
#> 
#> $nb_control$maxit_inner
#> [1] 100
#> 
#> $nb_control$theta_max
#> [1] 1e+06
#> 
#> $nb_control$trace
#> [1] FALSE
#> 
#> 
#> $maxit_logit
#> [1] 100
#> 
#> $epsilon_logit
#> [1] 1e-12
#> 
#> $hessian
#> [1] "analytic"
#> 
#> $trace
#> [1] TRUE
#> 
```
