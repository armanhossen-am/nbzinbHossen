#' Methods for nbreg objects
#'
#' Standard S3 methods for fitted negative binomial regressions returned by
#' [nbreg()].
#'
#' @param object,x An object of class `"nbreg"` (or `"summary.nbreg"` for
#'   the print method of the summary).
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"deviance"` (default), `"pearson"` or
#'   `"response"`. For `predict()`: `"link"` (linear predictor
#'   \eqn{\eta}) or `"response"` (mean \eqn{\mu}).
#' @param newdata Optional data frame in which to look for covariates (and
#'   offset variables) for prediction. If omitted, the fitted values are used.
#' @param signif.stars Logical; show significance stars?
#' @param ... Further arguments (currently unused).
#'
#' @return `coef()` a named vector; `vcov()` a matrix; `logLik()` an object
#'   of class `"logLik"` with `df` = number of coefficients + 1 (for
#'   \eqn{\theta}), so that [stats::AIC()] and [stats::BIC()] work;
#'   `fitted()`, `residuals()` and `predict()` numeric vectors; `nobs()` the
#'   number of observations with positive weight; `summary()` an object of
#'   class `"summary.nbreg"` with a coefficient table of estimates, standard
#'   errors, z values and p-values.
#'
#' @name nbreg-methods
#' @examples
#' data("bioChemists", package = "pscl")
#' fit <- nbreg(art ~ fem + kid5 + ment, data = bioChemists)
#' fit
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#' logLik(fit)
#' BIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(fem = "Women", kid5 = 0, ment = 10),
#'         type = "response")
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Negative binomial (NB2) regression, log link, fitted by IRLS\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  (SE", format(x$SE.theta, digits = digits), ")\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L),
      "on", length(x$coefficients) + 1L, "df\n")
  cat("Iterations:", x$iter, if (!x$converged) "(NOT converged)", "\n\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  zval <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = zval,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(zval)))
  out <- list(call = object$call, coefficients = tab, theta = object$theta,
              SE.theta = object$SE.theta, loglik = object$loglik,
              df = length(est) + 1L, aic = stats::AIC(object),
              iter = object$iter, converged = object$converged,
              deviance.resid = stats::residuals(object, type = "deviance"))
  class(out) <- "summary.nbreg"
  out
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Deviance residuals:\n")
  q <- stats::quantile(x$deviance.resid)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCoefficients (NB2, log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, na.print = "NA")
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Std. Err.:", format(x$SE.theta, digits = digits), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      x$df, "df,  AIC:", format(x$aic, digits = digits + 2L), "\n")
  cat("Iterations (IRLS / Newton-Raphson alternations):", x$iter,
      if (!x$converged) "(NOT converged)", "\n\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients) + 1L,
            nobs = stats::nobs(object), class = "logLik")
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) sum(object$weights > 0)

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("deviance", "pearson", "response"),
                            ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  th <- object$theta
  a <- object$weights
  switch(type,
         response = y - mu,
         # (y - mu) / sqrt(V(mu)) with V(mu) = mu + mu^2 / theta
         pearson = (y - mu) * sqrt(a) / sqrt(mu + mu^2 / th),
         deviance = {
           ylogy <- ifelse(y > 0, y * log(y / mu), 0)
           d <- 2 * a * (ylogy - (y + th) * log((y + th) / (mu + th)))
           sign(y - mu) * sqrt(pmax(d, 0))
         })
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    mt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(mt, newdata, na.action = stats::na.pass,
                             xlev = object$xlevels)
    X <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts)
    off <- stats::model.offset(mf)
    if (is.null(off)) off <- rep(0, nrow(X))
    if (!is.null(object$call$offset)) {
      off <- off + eval(object$call$offset, newdata,
                        environment(object$terms))
    }
    eta <- drop(X %*% object$coefficients) + off
  }
  if (type == "link") eta else exp(eta)
}
