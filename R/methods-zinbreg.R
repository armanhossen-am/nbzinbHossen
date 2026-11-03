#' Methods for zinbreg objects
#'
#' Standard S3 methods for fitted zero-inflated negative binomial
#' regressions returned by [zinbreg()].
#'
#' @param object,x An object of class `"zinbreg"` (or `"summary.zinbreg"`
#'   for the print method of the summary).
#' @param model For `coef()` and `vcov()`: `"full"` (count and zero
#'   coefficients, named `count_*` and `zero_*` as in `pscl::zeroinfl()`),
#'   `"count"` or `"zero"`.
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"pearson"` (default) or `"response"`.
#'   For `predict()`: `"response"` (marginal mean \eqn{(1-\pi)\mu}),
#'   `"count"` (NB mean \eqn{\mu}), `"zero"` (structural-zero probability
#'   \eqn{\pi}) or `"prob"` (matrix of \eqn{P(Y = k)} for `k` in `at`).
#' @param newdata Optional data frame in which to look for covariates.
#' @param at For `type = "prob"`: the counts at which to evaluate the
#'   probabilities (default `0:max(y)`).
#' @param signif.stars Logical; show significance stars?
#' @param ... Further arguments (currently unused).
#'
#' @return `coef()` a named vector; `vcov()` a matrix; `logLik()` an object
#'   of class `"logLik"` with `df` = number of count and zero coefficients
#'   + 1 (for \eqn{\theta}); `fitted()`, `residuals()` and `predict()`
#'   numeric vectors (a matrix for `type = "prob"`); `nobs()` the number of
#'   observations with positive weight; `summary()` an object of class
#'   `"summary.zinbreg"` with coefficient tables for both model parts.
#'
#' @name zinbreg-methods
#' @examples
#' data("bioChemists", package = "pscl")
#' fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
#' coef(fit, model = "zero")
#' AIC(fit)
#' nd <- data.frame(fem = "Men", kid5 = 0, ment = c(0, 10, 40))
#' predict(fit, newdata = nd, type = "response")
#' predict(fit, newdata = nd, type = "zero")
#' # expected number of zeros
#' sum(predict(fit, type = "prob", at = 0))
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Count model coefficients (negbin with log link):\n")
  print.default(format(x$coefficients$count, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      attr(stats::logLik(x), "df"), "df\n")
  cat("EM iterations:", x$iter, if (!x$converged) "(NOT converged)", "\n\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  se <- sqrt(diag(object$vcov_full))
  ctab <- function(est, se) {
    zval <- est / se
    cbind(Estimate = est, `Std. Error` = se, `z value` = zval,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(zval)))
  }
  count <- ctab(c(object$coefficients$count, `Log(theta)` = log(object$theta)),
                se[c(seq_len(p), p + q + 1L)])
  zero <- ctab(object$coefficients$zero, se[p + seq_len(q)])
  rownames(count) <- c(names(object$coefficients$count), "Log(theta)")
  rownames(zero) <- names(object$coefficients$zero)
  out <- list(call = object$call, coefficients = list(count = count, zero = zero),
              theta = object$theta, SE.theta = object$SE.theta,
              SE.logtheta = object$SE.logtheta, loglik = object$loglik,
              df = attr(stats::logLik(object), "df"), aic = stats::AIC(object),
              iter = object$iter, converged = object$converged,
              pearson.resid = stats::residuals(object, type = "pearson"))
  class(out) <- "summary.zinbreg"
  out
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson.resid)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCount model coefficients (negbin with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits,
                      signif.stars = signif.stars, na.print = "NA")
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits,
                      signif.stars = signif.stars, na.print = "NA")
  cat("\nTheta =", format(x$theta, digits = digits),
      "  (SE", format(x$SE.theta, digits = digits), ")\n")
  cat("Number of EM iterations:", x$iter,
      if (!x$converged) "(NOT converged)", "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      x$df, "df,  AIC:", format(x$aic, digits = digits + 2L), "\n\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         count = cf$count,
         zero = cf$zero,
         full = c(stats::setNames(cf$count, paste0("count_", names(cf$count))),
                  stats::setNames(cf$zero, paste0("zero_", names(cf$zero)))))
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  vc <- object$vcov
  switch(model,
         full = vc,
         count = {
           v <- vc[seq_len(p), seq_len(p), drop = FALSE]
           dimnames(v) <- list(names(object$coefficients$count),
                               names(object$coefficients$count))
           v
         },
         zero = {
           v <- vc[p + seq_len(q), p + seq_len(q), drop = FALSE]
           dimnames(v) <- list(names(object$coefficients$zero),
                               names(object$coefficients$zero))
           v
         })
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  df <- length(object$coefficients$count) + length(object$coefficients$zero) + 1L
  structure(object$loglik, df = df, nobs = stats::nobs(object),
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) sum(object$weights > 0)

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  res <- object$y - object$fitted.values
  if (type == "response") return(res)
  # Var(Y) = (1 - pi) mu (1 + mu (pi + 1/theta))
  mu <- object$mu
  pi <- object$pi
  v <- (1 - pi) * mu * (1 + mu * (pi + 1 / object$theta))
  res / sqrt(v)
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero", "prob"),
                            at = NULL, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    mu <- object$mu
    pi <- object$pi
  } else {
    mf <- stats::model.frame(stats::delete.response(object$terms$full),
                             newdata, na.action = stats::na.pass,
                             xlev = object$levels)
    mtX <- stats::delete.response(object$terms$count)
    X <- stats::model.matrix(mtX, mf, contrasts.arg = object$contrasts$count)
    Z <- stats::model.matrix(object$terms$zero, mf,
                             contrasts.arg = object$contrasts$zero)
    off <- count_offset(mtX, mf)
    if (!is.null(object$call$offset)) {
      off <- off + eval(object$call$offset, newdata,
                        environment(object$terms$full))
    }
    mu <- exp(drop(X %*% object$coefficients$count) + off)
    pi <- stats::plogis(drop(Z %*% object$coefficients$zero))
  }
  switch(type,
         response = (1 - pi) * mu,
         count = mu,
         zero = pi,
         prob = {
           if (is.null(at)) at <- 0:max(object$y)
           out <- vapply(at, function(k) {
             (1 - pi) * stats::dnbinom(k, size = object$theta, mu = mu) +
               pi * (k == 0)
           }, numeric(length(mu)))
           out <- matrix(out, nrow = length(mu))
           colnames(out) <- at
           out
         })
}
