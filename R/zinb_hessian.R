# Analytic Hessian of the observed-data ZINB log-likelihood with respect to
# (beta, gamma, phi), phi = log(theta).
#
# Notation: eta_i = x_i'beta + o_i, mu_i = exp(eta_i), zeta_i = z_i'gamma,
# pi_i = plogis(zeta_i), th = theta. Each l_i depends on beta only through
# eta_i and on gamma only through zeta_i, so with per-observation second
# derivatives l_ee, l_ez, l_ep, l_zz, l_zp, l_pp the blocks are
#   H_bb = X' diag(a l_ee) X,  H_bg = X' diag(a l_ez) Z,  H_bp = X'(a l_ep),
#   H_gg = Z' diag(a l_zz) Z,  H_gp = Z'(a l_zp),       H_pp = sum(a l_pp).
#
# y_i > 0:  l_i = log(1 - pi_i) + log f_NB(y_i; mu_i, th)
#   l_ee = -th mu (th + y) / (th + mu)^2
#   l_ep =  th mu (y - mu) / (th + mu)^2
#   l_pp =  th^2 H_i + th U_i       (U_i, H_i as in Section 3.1(b))
#   l_zz = -pi (1 - pi),  l_ez = l_zp = 0
#
# y_i = 0:  l_i = log(L_i),  L_i = pi_i + (1 - pi_i) f0_i,  f0_i = exp(g_i),
#   g = th log{th / (th + mu)}, with derivatives
#   g_e  = -th mu / (th + mu)
#   g_p  =  th [log{th / (th + mu)} + mu / (th + mu)]
#   g_ee = -th^2 mu / (th + mu)^2
#   g_ep = -th mu^2 / (th + mu)^2
#   g_pp =  th mu^2 / (th + mu)^2 + g_p
#   With r = (1 - pi) f0 / L (= 1 - tau, and dr/dg = r (1 - r)):
#   l_ee = r g_ee + r (1 - r) g_e^2
#   l_ep = r g_ep + r (1 - r) g_e g_p
#   l_pp = r g_pp + r (1 - r) g_p^2
#   l_ez = -pi (1 - pi) f0 g_e / L^2      (since dr/dzeta = -pi (1 - pi) f0 / L^2)
#   l_zp = -pi (1 - pi) f0 g_p / L^2
#   l_zz = (1 - f0) pi (1 - pi) [(1 - 2 pi) L - pi (1 - pi) (1 - f0)] / L^2
zinb_hessian <- function(X, Z, y, offset, a, beta, gamma, theta) {
  th <- theta
  mu <- exp(drop(X %*% beta) + offset)
  pi <- stats::plogis(drop(Z %*% gamma))
  pq <- pi * (1 - pi)
  tm <- th + mu
  pos <- y > 0

  # y > 0 contributions
  U_i <- digamma(y + th) - digamma(th) + log(th) + 1 - log(tm) - (y + th) / tm
  H_i <- trigamma(y + th) - trigamma(th) + 1 / th - 2 / tm + (y + th) / tm^2
  l_ee <- -th * mu * (th + y) / tm^2
  l_ep <- th * mu * (y - mu) / tm^2
  l_pp <- th^2 * H_i + th * U_i
  l_zz <- -pq
  l_ez <- numeric(length(y))
  l_zp <- numeric(length(y))

  # y = 0 contributions
  z0 <- !pos
  if (any(z0)) {
    m0 <- mu[z0]
    t0 <- tm[z0]
    p0 <- pi[z0]
    pq0 <- pq[z0]
    g <- th * log(th / t0)
    f0 <- exp(g)
    L <- p0 + (1 - p0) * f0
    r <- (1 - p0) * f0 / L
    g_e <- -th * m0 / t0
    g_p <- th * (log(th / t0) + m0 / t0)
    g_ee <- -th^2 * m0 / t0^2
    g_ep <- -th * m0^2 / t0^2
    g_pp <- th * m0^2 / t0^2 + g_p
    l_ee[z0] <- r * g_ee + r * (1 - r) * g_e^2
    l_ep[z0] <- r * g_ep + r * (1 - r) * g_e * g_p
    l_pp[z0] <- r * g_pp + r * (1 - r) * g_p^2
    l_ez[z0] <- -pq0 * f0 * g_e / L^2
    l_zp[z0] <- -pq0 * f0 * g_p / L^2
    l_zz[z0] <- (1 - f0) * pq0 * ((1 - 2 * p0) * L - pq0 * (1 - f0)) / L^2
  }

  H_bb <- crossprod(X, X * (a * l_ee))
  H_bg <- crossprod(X, Z * (a * l_ez))
  H_bp <- crossprod(X, a * l_ep)
  H_gg <- crossprod(Z, Z * (a * l_zz))
  H_gp <- crossprod(Z, a * l_zp)
  H_pp <- sum(a * l_pp)
  rbind(cbind(H_bb, H_bg, H_bp),
        cbind(t(H_bg), H_gg, H_gp),
        cbind(t(H_bp), t(H_gp), H_pp))
}
