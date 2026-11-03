#' Control parameters for nbreg()
#'
#' @param epsilon Convergence tolerance for the outer loop: iteration stops
#'   when the log-likelihood changes by less than `epsilon` between two
#'   alternations of the beta and theta updates.
#' @param maxit Maximum number of outer alternations.
#' @param epsilon_inner Relative tolerance for the inner IRLS loop (beta) and
#'   the absolute tolerance on `log(theta)` for the Newton-Raphson loop.
#' @param maxit_inner Maximum number of iterations of each inner loop.
#' @param theta_max Upper bound for theta. If the estimate reaches it the
#'   data show no overdispersion (Poisson limit) and a warning is given.
#' @param trace Logical; print the log-likelihood and theta at each outer
#'   iteration?
#'
#' @return A list of control parameters.
#' @export
#' @examples
#' nbreg_control(epsilon = 1e-10, trace = TRUE)
nbreg_control <- function(epsilon = 1e-8, maxit = 100, epsilon_inner = 1e-10,
                          maxit_inner = 100, theta_max = 1e6, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0) stop("'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit < 1) stop("'maxit' must be >= 1")
  if (!is.numeric(epsilon_inner) || epsilon_inner <= 0) {
    stop("'epsilon_inner' must be > 0")
  }
  if (!is.numeric(maxit_inner) || maxit_inner < 1) {
    stop("'maxit_inner' must be >= 1")
  }
  if (!is.numeric(theta_max) || theta_max <= 0) stop("'theta_max' must be > 0")
  list(epsilon = epsilon, maxit = as.integer(maxit),
       epsilon_inner = epsilon_inner, maxit_inner = as.integer(maxit_inner),
       theta_max = theta_max, trace = isTRUE(trace))
}

#' Negative binomial regression by IRLS
#'
#' Fits the NB2 regression model \eqn{\log \mu_i = x_i^\top \beta + o_i},
#' \eqn{\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta}, by alternating
#' (a) iteratively reweighted least squares for \eqn{\beta} with \eqn{\theta}
#' fixed and (b) Newton-Raphson for \eqn{\log\theta} with \eqn{\beta} fixed,
#' until the log-likelihood changes by less than `control$epsilon`.
#' No general-purpose optimizer is used. The results agree with
#' `MASS::glm.nb()`.
#'
#' @param formula A model formula, e.g. `y ~ x1 + x2`. Offsets can be given
#'   in the formula as `offset(log(exposure))`.
#' @param data A data frame containing the variables in the model.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values (default from
#'   `getOption("na.action")`).
#' @param weights Optional non-negative prior weights \eqn{a_i}.
#' @param offset Optional offset vector (added to the linear predictor).
#' @param control A list from [nbreg_control()].
#' @param model Logical; keep the model frame in the result?
#' @param ... Arguments passed to [nbreg_control()] when `control` is not
#'   given.
#'
#' @return An object of class `"nbreg"`, a list with components including
#'   `coefficients`, `vcov` (the covariance of \eqn{\hat\beta},
#'   \eqn{(X^\top \hat W X)^{-1}}, treating \eqn{\theta} as known),
#'   `theta`, `SE.theta`, `loglik`, `iter`, `converged`, `fitted.values`,
#'   `linear.predictors`, `y`, `x`, `weights`, `offset`, `df.residual` and
#'   `call`.
#'
#' @seealso [zinbreg()], [summary.nbreg()], [predict.nbreg()]
#' @export
#' @examples
#' data("bioChemists", package = "pscl")
#' fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#' summary(fit)
#' fit$theta
#' AIC(fit)
nbreg <- function(formula, data, subset, na.action, weights, offset,
                  control = nbreg_control(...), model = TRUE, ...) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  check_counts(y)
  X <- stats::model.matrix(mt, mf)
  n <- length(y)
  a <- stats::model.weights(mf)
  if (is.null(a)) a <- rep(1, n)
  if (any(a < 0)) stop("negative weights are not allowed", call. = FALSE)
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, n)
  if (length(off) != n) stop("'offset' has the wrong length", call. = FALSE)

  fit <- nb_fit(X, y, off, a, control)
  warn_nb_fit(fit, control)

  fit$y <- y
  fit$x <- X
  fit$weights <- a
  fit$offset <- off
  fit$df.residual <- sum(a > 0) - ncol(X)
  fit$call <- cl
  fit$terms <- mt
  fit$xlevels <- stats::.getXlevels(mt, mf)
  fit$contrasts <- attr(X, "contrasts")
  fit$control <- control
  if (model) fit$model <- mf
  class(fit) <- "nbreg"
  fit
}

# Matrix-level NB fit, Section 3.1. Also used as the count-part M-step of
# zinbreg() with prior weights a_i = 1 - tau_i, so it accepts starting values.
#
# (c) Outer loop: start from a Poisson IRLS fit and the moment estimate
#     theta_0 = sum(a) / sum(a (y/mu - 1)^2), then alternate
#     (a) IRLS for beta and (b) Newton-Raphson for log(theta)
#     until |l_new - l_old| < epsilon.
# (d) At convergence: Var(beta) = (X'WX)^{-1}, SE(theta) = {-H(theta)}^{-1/2}.
nb_fit <- function(X, y, offset, a, control, beta = NULL, theta = NULL) {
  if (is.null(beta)) beta <- poisson_irls(X, y, offset, a)
  mu <- exp(drop(X %*% beta) + offset)
  if (is.null(theta)) {
    theta <- sum(a) / sum(a * (y / mu - 1)^2)
    if (!is.finite(theta) || theta <= 0) theta <- 1
    theta <- min(theta, control$theta_max)
  }
  ll <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  at_max <- FALSE
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    fb <- nb_irls_beta(X, y, offset, a, theta, beta,
                       maxit = control$maxit_inner, tol = control$epsilon_inner)
    beta <- fb$beta
    mu <- fb$mu
    ft <- nb_theta_newton(y, mu, a, theta, maxit = control$maxit_inner,
                          tol = control$epsilon_inner,
                          theta_max = control$theta_max)
    theta <- ft$theta
    at_max <- ft$at_max
    ll_new <- ft$loglik
    if (control$trace) {
      message(sprintf("nbreg iter %3d: logLik = %.10f, theta = %.6f",
                      iter, ll_new, theta))
    }
    dll <- ll_new - ll
    ll <- ll_new
    if (abs(dll) < control$epsilon) {
      converged <- TRUE
      break
    }
  }

  eta <- drop(X %*% beta) + offset
  mu <- exp(eta)
  names(beta) <- colnames(X)
  w <- a * theta * mu / (theta + mu)          # IRLS weights at convergence
  H <- nb_theta_hess(y, mu, theta, a)
  list(coefficients = beta,
       vcov = xtwx_inverse(X, w),
       theta = theta,
       SE.theta = if (H < 0) 1 / sqrt(-H) else NA_real_,
       loglik = ll,
       iter = iter,
       converged = converged,
       theta_at_max = at_max,
       fitted.values = mu,
       linear.predictors = eta,
       irls_weights = w)
}

warn_nb_fit <- function(fit, control) {
  if (fit$theta_at_max) {
    warning("theta reached 'theta_max' (", format(control$theta_max),
            "): no evidence of overdispersion, a Poisson model may suffice",
            call. = FALSE)
  } else if (fit$theta > 1e4) {
    warning("theta estimate is very large (", format(fit$theta, digits = 4),
            "); the fit is close to the Poisson limit", call. = FALSE)
  }
  if (!fit$converged) {
    warning("nbreg did not converge in ", control$maxit, " iterations",
            call. = FALSE)
  }
  invisible(NULL)
}

check_counts <- function(y) {
  if (is.null(y)) stop("the formula must have a response", call. = FALSE)
  if (!is.numeric(y)) stop("the response must be numeric", call. = FALSE)
  if (any(y < 0)) stop("the response must be non-negative counts", call. = FALSE)
  if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
    stop("the response must be integer counts", call. = FALSE)
  }
  invisible(NULL)
}
