# Control parameters for nbreg()

Control parameters for nbreg()

## Usage

``` r
nbreg_control(
  epsilon = 1e-08,
  maxit = 100,
  epsilon_inner = 1e-10,
  maxit_inner = 100,
  theta_max = 1e+06,
  trace = FALSE
)
```

## Arguments

- epsilon:

  Convergence tolerance for the outer loop: iteration stops when the
  log-likelihood changes by less than `epsilon` between two alternations
  of the beta and theta updates.

- maxit:

  Maximum number of outer alternations.

- epsilon_inner:

  Relative tolerance for the inner IRLS loop (beta) and the absolute
  tolerance on `log(theta)` for the Newton-Raphson loop.

- maxit_inner:

  Maximum number of iterations of each inner loop.

- theta_max:

  Upper bound for theta. If the estimate reaches it the data show no
  overdispersion (Poisson limit) and a warning is given.

- trace:

  Logical; print the log-likelihood and theta at each outer iteration?

## Value

A list of control parameters.

## Examples

``` r
nbreg_control(epsilon = 1e-10, trace = TRUE)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
#> [1] 100
#> 
#> $epsilon_inner
#> [1] 1e-10
#> 
#> $maxit_inner
#> [1] 100
#> 
#> $theta_max
#> [1] 1e+06
#> 
#> $trace
#> [1] TRUE
#> 
```
