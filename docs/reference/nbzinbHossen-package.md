# nbzinbHossen: NB and ZINB regression by IRLS and EM

Negative binomial regression is fitted with
[`nbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/nbreg.md)
by iteratively reweighted least squares alternated with Newton-Raphson
for the dispersion parameter. Zero-inflated negative binomial regression
is fitted with
[`zinbreg()`](https://armanhossen-am.github.io/nbzinbHossen/reference/zinbreg.md)
by an EM algorithm whose M-steps are IRLS fits. No general-purpose
optimizer (such as [`optim()`](https://rdrr.io/r/stats/optim.html)) is
used.

## See also

Useful links:

- <https://github.com/armanhossen-am/nbzinbHossen>

- <https://armanhossen-am.github.io/nbzinbHossen/>

- Report bugs at <https://github.com/armanhossen-am/nbzinbHossen/issues>

## Author

**Maintainer**: Arman Hossen <hossena@miamioh.edu>

Authors:

- Arman Hossen <hossena@miamioh.edu>
