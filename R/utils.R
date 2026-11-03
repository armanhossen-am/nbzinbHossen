# Internal numerical routines shared by nbreg() and zinbreg().
#
# Nothing in this file calls optim(), glm(), lm() or any other fitting or
# optimization routine: every update is an explicit weighted least-squares
# solve (via qr()) or an explicit Newton-Raphson step.

# Weighted least squares: returns argmin_b sum_i w_i (z_i - x_i' b)^2,
# i.e. b = (X'WX)^{-1} X'Wz, computed stably from the QR decomposition of
# W^{1/2} X instead of forming X'WX explicitly.
wls_solve <- function(X, w, z) {
  sw <- sqrt(w)
  qrx <- qr(X * sw)
  if (qrx$rank < ncol(X)) {
    stop("design matrix is rank deficient; remove collinear covariates",
         call. = FALSE)
  }
  drop(qr.coef(qrx, z * sw))
}

# Inverse of X'WX via its Cholesky factor (used for the NB covariance).
xtwx_inverse <- function(X, w) {
  out <- chol2inv(chol(crossprod(X * sqrt(w))))
  dimnames(out) <- list(colnames(X), colnames(X))
  out
}

# NB2 log-likelihood contributions, Section 2.1:
#   log Gamma(y + theta) - log Gamma(theta) - log y! + theta log theta
#   + y log mu - (theta + y) log(theta + mu).
# dnbinom(log = TRUE) evaluates exactly this expression, with care for
# large theta and y = 0.
nb_loglik_i <- function(y, mu, theta) {
  stats::dnbinom(y, size = theta, mu = mu, log = TRUE)
}

nb_loglik <- function(y, mu, theta, a) {
  sum(a * nb_loglik_i(y, mu, theta))
}

# Poisson IRLS, used only to get starting values for beta.
# Log link, V(mu) = mu:  w_i = a_i mu_i,  z_i = (eta_i - o_i) + (y_i - mu_i)/mu_i.
poisson_irls <- function(X, y, offset, a, maxit = 50, tol = 1e-8) {
  mu <- y + 0.1                       # same starting values as glm(family = poisson)
  eta <- log(mu)
  beta <- rep(0, ncol(X))
  for (it in seq_len(maxit)) {
    w <- a * mu
    z <- (eta - offset) + (y - mu) / mu
    beta_new <- wls_solve(X, w, z)
    eta <- drop(X %*% beta_new) + offset
    mu <- exp(eta)
    done <- max(abs(beta_new - beta)) < tol * (max(abs(beta)) + tol)
    beta <- beta_new
    if (done) break
  }
  beta
}

# (a) IRLS for beta with theta held fixed, Section 3.1(a).
# Log link, V(mu) = mu + mu^2/theta, so (d mu/d eta)^2 = mu^2 and
#   w_i = a_i mu_i^2 / V(mu_i) = a_i theta mu_i / (theta + mu_i),
#   z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i,
#   beta <- (X'WX)^{-1} X'Wz.
# The log link is not canonical for the NB, so Fisher scoring is not
# guaranteed to increase the likelihood; a step-halving safeguard is used.
nb_irls_beta <- function(X, y, offset, a, theta, beta, maxit, tol) {
  eta <- drop(X %*% beta) + offset
  mu <- exp(eta)
  ll <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    w <- a * theta * mu / (theta + mu)
    z <- (eta - offset) + (y - mu) / mu
    beta_new <- wls_solve(X, w, z)
    step <- beta_new - beta
    for (half in 0:30) {
      eta_new <- drop(X %*% (beta + step)) + offset
      mu_new <- exp(eta_new)
      ll_new <- nb_loglik(y, mu_new, theta, a)
      if (is.finite(ll_new) && ll_new >= ll - 1e-10 * abs(ll)) break
      step <- step / 2
    }
    beta <- beta + step
    eta <- eta_new
    mu <- mu_new
    dll <- ll_new - ll
    ll <- ll_new
    if (abs(dll) < tol * (abs(ll) + tol)) {
      converged <- TRUE
      break
    }
  }
  list(beta = beta, eta = eta, mu = mu, loglik = ll, iter = it,
       converged = converged)
}

# Score U(theta) and second derivative H(theta) of the NB log-likelihood in
# theta with mu held fixed, Section 3.1(b):
#   U = sum a_i [psi(y_i + th) - psi(th) + log th + 1 - log(th + mu_i)
#                - (y_i + th)/(th + mu_i)]
#   H = sum a_i [psi'(y_i + th) - psi'(th) + 1/th - 2/(th + mu_i)
#                + (y_i + th)/(th + mu_i)^2]
nb_theta_score <- function(y, mu, theta, a) {
  sum(a * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
             log(theta + mu) - (y + theta) / (theta + mu)))
}

nb_theta_hess <- function(y, mu, theta, a) {
  sum(a * (trigamma(y + theta) - trigamma(theta) + 1 / theta -
             2 / (theta + mu) + (y + theta) / (theta + mu)^2))
}

# (b) Newton-Raphson for theta with beta (hence mu) held fixed, Section 3.1(b).
# Works on phi = log(theta) so that theta stays positive:
#   U_phi = theta U,  H_phi = theta^2 H + theta U,  phi <- phi - U_phi / H_phi.
# If H_phi >= 0 (log-likelihood not locally concave in phi) a unit step in
# the direction of U_phi is used instead; in both cases the step is halved
# until the log-likelihood does not decrease.
nb_theta_newton <- function(y, mu, a, theta, maxit, tol, theta_max) {
  phi <- log(theta)
  ll <- nb_loglik(y, mu, theta, a)
  at_max <- FALSE
  converged <- FALSE
  for (it in seq_len(maxit)) {
    th <- exp(phi)
    U <- nb_theta_score(y, mu, th, a)
    H <- nb_theta_hess(y, mu, th, a)
    U_phi <- th * U
    H_phi <- th^2 * H + th * U
    step <- if (H_phi < 0) -U_phi / H_phi else sign(U_phi)
    step <- max(min(step, 5), -5)     # never move theta by more than e^5 at once
    for (half in 0:30) {
      phi_new <- min(phi + step, log(theta_max))
      ll_new <- nb_loglik(y, mu, exp(phi_new), a)
      if (is.finite(ll_new) && ll_new >= ll - 1e-12 * abs(ll)) break
      step <- step / 2
    }
    moved <- abs(phi_new - phi)
    phi <- phi_new
    ll <- ll_new
    if (phi >= log(theta_max)) {
      at_max <- TRUE
      break
    }
    if (moved < tol) {
      converged <- TRUE
      break
    }
  }
  list(theta = exp(phi), loglik = ll, iter = it, converged = converged,
       at_max = at_max)
}

# IRLS for a logistic regression with a (possibly fractional) response
# tau_i in [0, 1] and prior weights a_i, used for the ZINB zero-part M-step,
# Section 3.2:
#   maximize sum a_i [tau_i log pi_i + (1 - tau_i) log(1 - pi_i)],
#   w_i = a_i pi_i (1 - pi_i),  z_i = z_i'gamma + (tau_i - pi_i) / {pi_i (1 - pi_i)}.
logit_irls <- function(Z, tau, a, gamma, maxit = 100, tol = 1e-10) {
  obj <- function(eta) {
    sum(a * (tau * stats::plogis(eta, log.p = TRUE) +
               (1 - tau) * stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)))
  }
  eta <- drop(Z %*% gamma)
  ll <- obj(eta)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    p <- stats::plogis(eta)
    v <- pmax(p * (1 - p), 1e-12)     # guard against p = 0 or 1 exactly
    w <- a * v
    z <- eta + (tau - p) / v
    gamma_new <- wls_solve(Z, w, z)
    step <- gamma_new - gamma
    for (half in 0:30) {
      eta_new <- drop(Z %*% (gamma + step))
      ll_new <- obj(eta_new)
      if (is.finite(ll_new) && ll_new >= ll - 1e-10 * abs(ll)) break
      step <- step / 2
    }
    gamma <- gamma + step
    eta <- eta_new
    dll <- ll_new - ll
    ll <- ll_new
    if (abs(dll) < tol * (abs(ll) + tol)) {
      converged <- TRUE
      break
    }
  }
  list(gamma = gamma, loglik = ll, iter = it, converged = converged)
}
