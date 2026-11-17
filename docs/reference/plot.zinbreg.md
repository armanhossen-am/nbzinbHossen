# Plot the EM convergence of a zinbreg fit

Shows the observed-data log-likelihood at every EM iteration of
[`zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md).
EM never decreases the log-likelihood, so the curve climbs monotonically
to its final value.

## Usage

``` r
# S3 method for class 'zinbreg'
plot(x, which = c("loglik", "gap"), reference = NULL, ...)
```

## Arguments

- x:

  An object of class `"zinbreg"`.

- which:

  `"loglik"` (default) plots the log-likelihood against the iteration.
  `"gap"` plots the distance to the final value, \\(\ell\_{final} -
  \ell_t\\), on a log scale; a straight line there means linear
  convergence.

- reference:

  Optional log-likelihood of another fit of the same model, for example
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) of
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html),
  drawn as a dashed line when `which = "loglik"`.

- ...:

  Further graphical arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main` or `col`.

## Value

Invisibly, a data frame with columns `iteration` and `loglik`. Iteration
0 is the log-likelihood at the starting values.

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
plot(fit)

plot(fit, which = "gap")
```
