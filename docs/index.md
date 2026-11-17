# nbzinbHossen

**nbzinbHossen** fits

- **negative binomial (NB2) regression** with
  [`nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md)
  (an alternative to
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)), using
  iteratively reweighted least squares (IRLS) for the coefficients
  alternated with Newton–Raphson on log(θ); and
- **zero-inflated negative binomial (ZINB) regression** with
  [`zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md)
  (an alternative to `pscl::zeroinfl(..., dist = "negbin")`), using an
  EM algorithm whose M-steps are IRLS fits (a fractional-response
  logistic regression for the zero part and a weighted NB regression for
  the count part).

No general-purpose optimizer
([`optim()`](https://rdrr.io/r/stats/optim.html),
[`nlm()`](https://rdrr.io/r/stats/nlm.html), …) and no
[`glm()`](https://rdrr.io/r/stats/glm.html)/[`lm()`](https://rdrr.io/r/stats/lm.html)
fitting code is used. ZINB standard errors come from an **analytic**
Hessian of the observed-data log-likelihood.

Website (function reference + project report):
<https://armanhossen-am.github.io/nbzinbHossen/>

## Installation

\
`# install.packages("devtools")`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"armanhossen-am/nbzinbHossen"``)`

## Usage

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`nbzinbHossen`](https://github.com/armanhossen-am/nbzinbHossen)`)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"bioChemists"``, package ``=`` ``"pscl"``)`\
\
`# Negative binomial regression`\
`nb`` ``<-`` `[`nbreg`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md)`(``art`` ``~`` ``fem`` ``+`` ``mar`` ``+`` ``kid5`` ``+`` ``phd`` ``+`` ``ment``, data ``=`` ``bioChemists``)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``nb``)`\
`nb``$``theta``                     ``# 2.264, same as MASS::glm.nb()`\
\
`# Zero-inflated NB: count covariates | zero-inflation covariates`\
`zinb`` ``<-`` `[`zinbreg`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md)`(``art`` ``~`` ``fem`` ``+`` ``mar`` ``+`` ``kid5`` ``+`` ``phd`` ``+`` ``ment`` ``|`` ``ment``,`\
`                data ``=`` ``bioChemists``)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``zinb``)`\
[`AIC`](https://rdrr.io/r/stats/AIC.html)`(``nb``, ``zinb``)`\
\
`# Predictions: marginal mean, NB mean, and structural-zero probability`\
`nd`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``fem ``=`` ``"Women"``, mar ``=`` ``"Married"``, kid5 ``=`` ``1``, phd ``=`` ``3``, ment ``=`` ``10``)`\
[`predict`](https://rdrr.io/r/stats/predict.html)`(``zinb``, newdata ``=`` ``nd``, type ``=`` ``"response"``)`\
[`predict`](https://rdrr.io/r/stats/predict.html)`(``zinb``, newdata ``=`` ``nd``, type ``=`` ``"count"``)`\
[`predict`](https://rdrr.io/r/stats/predict.html)`(``zinb``, newdata ``=`` ``nd``, type ``=`` ``"zero"``)`

Both fitting functions return S3 objects with `print`, `summary`,
`coef`, `vcov`, `logLik` (so
[`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
work), `fitted`, `residuals`, `nobs` and `predict` methods. Use
[`nbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg_control.md)
/
[`zinbreg_control()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg_control.md)
to change tolerances, iteration limits, or to print the iteration
history (`trace = TRUE`).

## Agreement with the reference implementations

On [`pscl::bioChemists`](https://rdrr.io/pkg/pscl/man/bioChemists.html)
and [`AER::NMES1988`](https://rdrr.io/pkg/AER/man/NMES1988.html) the
maximum absolute differences are about 1e-7 from
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) for NB, and
below 3e-5 from
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) for
ZINB. The tolerances required by the project are 1e-5 and 1e-3. The full
tables are in `inst/benchmark/` and in the report.

## R CMD check

[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
gives 0 errors, 0 warnings and 0 notes.
