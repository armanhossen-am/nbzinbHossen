# Simulation study (project Section 5.1): does zinbreg() recover known
# parameters? R = 200 data sets of size n = 500 from a ZINB model with two
# covariates in the count part and one in the zero part.
#
# Run from the package root:  Rscript inst/benchmark/simulation.R
# Writes inst/benchmark/simulation_results.rds.

library(nbzinbHossen)

run_simulation <- function(R = 200, n = 500, seed = 651) {
  beta <- c("(Intercept)" = 1.0, x1 = 0.5, x2 = -0.4)
  gamma <- c("(Intercept)" = -1.0, w = 0.8)
  theta <- 1.5
  truth <- c(paste0("count_", names(beta)), paste0("zero_", names(gamma)),
             "Log(theta)")
  truth <- stats::setNames(c(beta, gamma, log(theta)), truth)

  set.seed(seed)
  est <- se <- matrix(NA_real_, R, length(truth),
                      dimnames = list(NULL, names(truth)))
  conv <- logical(R)
  iters <- integer(R)
  for (r in seq_len(R)) {
    x1 <- stats::rnorm(n)
    x2 <- stats::rbinom(n, 1, 0.5)
    w <- stats::rnorm(n)
    mu <- exp(beta[1] + beta[2] * x1 + beta[3] * x2)
    pi <- stats::plogis(gamma[1] + gamma[2] * w)
    y <- ifelse(stats::runif(n) < pi, 0,
                stats::rnbinom(n, size = theta, mu = mu))
    dat <- data.frame(y, x1, x2, w)
    fit <- tryCatch(suppressWarnings(zinbreg(y ~ x1 + x2 | w, data = dat)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    est[r, ] <- c(coef(fit), log(fit$theta))
    se[r, ] <- sqrt(diag(fit$vcov_full))
    conv[r] <- fit$converged
    iters[r] <- fit$iter
  }

  ok <- stats::complete.cases(est) & conv
  est <- est[ok, , drop = FALSE]
  se <- se[ok, , drop = FALSE]
  lower <- est - stats::qnorm(0.975) * se
  upper <- est + stats::qnorm(0.975) * se
  cover <- sweep(lower, 2, truth, "<=") & sweep(upper, 2, truth, ">=")

  summary_table <- data.frame(
    parameter = names(truth),
    true = unname(truth),
    mean_est = colMeans(est),
    bias = colMeans(est) - truth,
    mc_sd = apply(est, 2, stats::sd),
    mean_se = colMeans(se),
    coverage = colMeans(cover),
    row.names = NULL
  )
  list(summary = summary_table, R = R, n = n, n_ok = sum(ok),
       median_iter = stats::median(iters[ok]), seed = seed)
}

if (sys.nframe() == 0L) {
  sim <- run_simulation()
  print(sim$summary, digits = 3)
  cat("Converged fits:", sim$n_ok, "of", sim$R,
      "; median EM iterations:", sim$median_iter, "\n")
  saveRDS(sim, file.path("inst", "benchmark", "simulation_results.rds"))
}
