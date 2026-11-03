#' nbzinbHossen: NB and ZINB regression by IRLS and EM
#'
#' Negative binomial regression is fitted with [nbreg()] by iteratively
#' reweighted least squares alternated with Newton-Raphson for the
#' dispersion parameter. Zero-inflated negative binomial regression is
#' fitted with [zinbreg()] by an EM algorithm whose M-steps are IRLS fits.
#' No general-purpose optimizer (such as `optim()`) is used.
#'
#' @keywords internal
"_PACKAGE"
