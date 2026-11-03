#' Control parameters for zinbreg()
#'
#' @param epsilon Convergence tolerance for the EM algorithm: iteration stops
#'   when the observed-data log-likelihood changes by less than `epsilon`.
#' @param maxit Maximum number of EM iterations.
#' @param nb_control Control list from [nbreg_control()] used for the initial
#'   NB fit and for the count-part M-step.
#' @param maxit_logit,epsilon_logit Maximum iterations and relative
#'   tolerance of the IRLS logistic regression in the zero-part M-step.
#' @param hessian How to compute the observed information for the standard
#'   errors: `"analytic"` (default, closed-form second derivatives) or
#'   `"numerical"` (`numDeriv::hessian()` of the observed-data
#'   log-likelihood; requires the numDeriv package).
#' @param trace Logical; print the observed-data log-likelihood at each EM
#'   iteration?
#'
#' @return A list of control parameters.
#' @export
#' @examples
#' zinbreg_control(epsilon = 1e-10, trace = TRUE)
zinbreg_control <- function(epsilon = 1e-10, maxit = 5000,
                            nb_control = nbreg_control(epsilon = 1e-10),
                            maxit_logit = 100, epsilon_logit = 1e-12,
                            hessian = c("analytic", "numerical"),
                            trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0) stop("'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit < 1) stop("'maxit' must be >= 1")
  hessian <- match.arg(hessian)
  list(epsilon = epsilon, maxit = as.integer(maxit), nb_control = nb_control,
       maxit_logit = as.integer(maxit_logit), epsilon_logit = epsilon_logit,
       hessian = hessian, trace = isTRUE(trace))
}

#' Zero-inflated negative binomial regression by EM
#'
#' Fits the ZINB model
#' \deqn{P(Y_i = 0) = \pi_i + (1 - \pi_i) f_{NB}(0; \mu_i, \theta), \quad
#'       P(Y_i = y) = (1 - \pi_i) f_{NB}(y; \mu_i, \theta),\ y \ge 1,}
#' with \eqn{\log\mu_i = x_i^\top\beta + o_i} and
#' \eqn{\mathrm{logit}\,\pi_i = z_i^\top\gamma}, by an EM algorithm. The
#' E-step computes the posterior probability \eqn{\tau_i} that observation
#' \eqn{i} is a structural zero; the M-step fits a logistic regression with
#' fractional response \eqn{\tau_i} by IRLS (zero part) and a weighted NB
#' regression with prior weights \eqn{1 - \tau_i} by [nbreg()]'s IRLS /
#' Newton-Raphson algorithm (count part). Standard errors come from the
#' observed information of the observed-data log-likelihood in
#' \eqn{(\beta, \gamma, \log\theta)}. The results agree with
#' `pscl::zeroinfl(..., dist = "negbin")`.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before `|`, zero-model covariates after it. Without `|` the
#'   same covariates are used in both parts. Offsets (in the formula or via
#'   `offset`) apply to the count part only.
#' @inheritParams nbreg
#' @param control A list from [zinbreg_control()].
#' @param ... Arguments passed to [zinbreg_control()] when `control` is not
#'   given.
#'
#' @return An object of class `"zinbreg"`, a list with components including
#'   `coefficients` (a list with `count` and `zero`), `vcov` (covariance of
#'   all count and zero coefficients), `vcov_full` (also including
#'   `log(theta)`), `theta`, `SE.logtheta`, `SE.theta`, `loglik`, `iter`,
#'   `converged`, `loglik_trace`, `fitted.values` (\eqn{(1-\pi_i)\mu_i}),
#'   `mu`, `pi`, `y`, `x`, `z`, `weights`, `offset`, `df.residual` and
#'   `call`.
#'
#' @seealso [nbreg()], [summary.zinbreg()], [predict.zinbreg()]
#' @export
#' @examples
#' data("bioChemists", package = "pscl")
#' fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
#'                data = bioChemists)
#' summary(fit)
#' head(predict(fit, type = "zero"))
zinbreg <- function(formula, data, subset, na.action, weights, offset,
                    control = zinbreg_control(...), model = TRUE, ...) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  # Two-part formula: y ~ count | zero. Without "|", reuse the count part.
  ff <- Formula::as.Formula(formula)
  if (length(ff)[2L] < 2L) {
    ff <- Formula::as.Formula(stats::formula(ff, rhs = 1L),
                              stats::formula(ff, lhs = 0L, rhs = 1L))
  } else if (length(ff)[2L] > 2L) {
    stop("formula must have at most two parts: y ~ count | zero",
         call. = FALSE)
  }
  mf$formula <- ff
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  mtX <- stats::terms(ff, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(ff, data = data, rhs = 2L))
  y <- stats::model.response(mf, "numeric")
  check_counts(y)
  if (all(y > 0)) {
    stop("the response has no zeros; a zero-inflated model is not identified",
         call. = FALSE)
  }
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  n <- length(y)
  a <- stats::model.weights(mf)
  if (is.null(a)) a <- rep(1, n)
  if (any(a < 0)) stop("negative weights are not allowed", call. = FALSE)
  off <- count_offset(mtX, mf)

  fit <- zinb_fit(X, Z, y, off, a, control)
  if (!fit$converged) {
    warning("zinbreg EM did not converge in ", control$maxit, " iterations",
            call. = FALSE)
  }
  if (fit$theta > 1e4) {
    warning("theta estimate is very large (", format(fit$theta, digits = 4),
            "); the count part is close to Poisson", call. = FALSE)
  }

  fit$y <- y
  fit$x <- X
  fit$z <- Z
  fit$weights <- a
  fit$offset <- off
  fit$df.residual <- sum(a > 0) - ncol(X) - ncol(Z) - 1L
  fit$call <- cl
  fit$formula <- ff
  fit$terms <- list(count = mtX, zero = mtZ, full = mt)
  fit$levels <- stats::.getXlevels(mt, mf)
  fit$contrasts <- list(count = attr(X, "contrasts"),
                        zero = attr(Z, "contrasts"))
  fit$control <- control
  if (model) fit$model <- mf
  class(fit) <- "zinbreg"
  fit
}

# Count-part offset: offset() terms in the count formula plus the 'offset'
# argument (stored by model.frame() as "(offset)").
count_offset <- function(mtX, mf) {
  off <- rep(0, nrow(mf))
  idx <- attr(mtX, "offset")
  if (length(idx)) {
    vars <- as.character(attr(mtX, "variables"))[-1L]
    for (i in idx) off <- off + mf[[vars[i]]]
  }
  if (!is.null(mf[["(offset)"]])) off <- off + mf[["(offset)"]]
  off
}

# Observed-data ZINB log-likelihood, Section 2.2:
#   y_i = 0: log[pi_i + (1 - pi_i) f_NB(0; mu_i, theta)]
#   y_i > 0: log(1 - pi_i) + log f_NB(y_i; mu_i, theta)
zinb_loglik_i <- function(y, mu, pi, theta) {
  f0 <- (theta / (theta + mu))^theta
  ifelse(y == 0,
         log(pi + (1 - pi) * f0),
         log1p(-pi) + nb_loglik_i(y, mu, theta))
}

zinb_loglik <- function(y, mu, pi, theta, a) {
  sum(a * zinb_loglik_i(y, mu, pi, theta))
}

# Log-likelihood as a function of the stacked vector (beta, gamma, log theta);
# used for the optional numDeriv Hessian and in the tests.
zinb_loglik_par <- function(par, X, Z, y, offset, a) {
  p <- ncol(X)
  q <- ncol(Z)
  mu <- exp(drop(X %*% par[seq_len(p)]) + offset)
  pi <- stats::plogis(drop(Z %*% par[p + seq_len(q)]))
  zinb_loglik(y, mu, pi, exp(par[p + q + 1L]), a)
}

# EM algorithm, Section 3.2.
zinb_fit <- function(X, Z, y, offset, a, control) {
  nb_ctrl <- control$nb_control
  is0 <- as.numeric(y == 0)

  # Starting values: NB fit for (beta, theta); IRLS logistic regression of
  # I(y = 0) on Z for gamma.
  nb0 <- nb_fit(X, y, offset, a, nb_ctrl)
  beta <- nb0$coefficients
  theta <- nb0$theta
  gamma <- logit_irls(Z, is0, a, rep(0, ncol(Z)),
                      maxit = control$maxit_logit,
                      tol = control$epsilon_logit)$gamma

  mu <- exp(drop(X %*% beta) + offset)
  pi <- stats::plogis(drop(Z %*% gamma))
  ll <- zinb_loglik(y, mu, pi, theta, a)
  ll_trace <- ll
  converged <- FALSE
  monotone <- TRUE
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    # E-step: tau_i = P(D_i = 1 | y_i)
    #   = pi_i / {pi_i + (1 - pi_i) f_NB(0; mu_i, theta)} if y_i = 0, else 0.
    f0 <- (theta / (theta + mu))^theta
    tau <- is0 * pi / (pi + (1 - pi) * f0)

    # M-step, zero part: logistic IRLS with fractional response tau_i.
    gamma <- logit_irls(Z, tau, a, gamma, maxit = control$maxit_logit,
                        tol = control$epsilon_logit)$gamma

    # M-step, count part: NB fit (Section 3.1) with prior weights a_i (1 - tau_i),
    # warm-started at the current (beta, theta).
    nbm <- nb_fit(X, y, offset, a * (1 - tau), nb_ctrl, beta = beta,
                  theta = theta)
    beta <- nbm$coefficients
    theta <- nbm$theta

    mu <- exp(drop(X %*% beta) + offset)
    pi <- stats::plogis(drop(Z %*% gamma))
    ll_new <- zinb_loglik(y, mu, pi, theta, a)
    ll_trace <- c(ll_trace, ll_new)
    # EM must never decrease the observed-data log-likelihood.
    if (ll_new < ll - 1e-8 * max(1, abs(ll))) monotone <- FALSE
    if (control$trace) {
      message(sprintf("zinbreg EM iter %4d: logLik = %.10f, theta = %.6f",
                      iter, ll_new, theta))
    }
    dll <- ll_new - ll
    ll <- ll_new
    if (abs(dll) < control$epsilon) {
      converged <- TRUE
      break
    }
  }
  if (!monotone) {
    warning("the EM log-likelihood decreased at some iteration; ",
            "check 'loglik_trace'", call. = FALSE)
  }

  names(beta) <- colnames(X)
  names(gamma) <- colnames(Z)
  p <- ncol(X)
  q <- ncol(Z)

  # Standard errors: invert the observed information of the observed-data
  # log-likelihood with respect to (beta, gamma, log theta).
  par <- c(beta, gamma, log(theta))
  H <- if (control$hessian == "analytic") {
    zinb_hessian(X, Z, y, offset, a, beta, gamma, theta)
  } else {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("hessian = \"numerical\" requires the numDeriv package",
           call. = FALSE)
    }
    numDeriv::hessian(zinb_loglik_par, par, X = X, Z = Z, y = y,
                      offset = offset, a = a)
  }
  nm <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
          "Log(theta)")
  vc <- tryCatch(solve(-H), error = function(e) {
    warning("observed information is singular; standard errors are NA",
            call. = FALSE)
    matrix(NA_real_, p + q + 1L, p + q + 1L)
  })
  dimnames(vc) <- list(nm, nm)
  se_logtheta <- sqrt(vc[p + q + 1L, p + q + 1L])

  list(coefficients = list(count = beta, zero = gamma),
       vcov = vc[seq_len(p + q), seq_len(p + q), drop = FALSE],
       vcov_full = vc,
       theta = theta,
       SE.logtheta = se_logtheta,
       SE.theta = theta * se_logtheta,          # delta method
       loglik = ll,
       loglik_trace = ll_trace,
       monotone = monotone,
       iter = iter,
       converged = converged,
       fitted.values = (1 - pi) * mu,
       mu = mu,
       pi = pi,
       tau = is0 * pi / (pi + (1 - pi) * (theta / (theta + mu))^theta))
}
