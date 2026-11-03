# Benchmark nbreg() / zinbreg() against MASS::glm.nb() / pscl::zeroinfl()
# on pscl::bioChemists and AER::NMES1988 (project Section 5).
#
# Run from the package root:  Rscript inst/benchmark/benchmark.R
# Writes inst/benchmark/benchmark_results.rds.

library(nbzinbHossen)
library(MASS)
library(pscl)
data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

datasets <- list(
  bioChemists = list(
    data = bioChemists,
    nb = art ~ fem + mar + kid5 + phd + ment,
    zinb = art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ),
  NMES1988 = list(
    data = NMES1988,
    nb = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance
  )
)

# Median elapsed time over 'reps' runs.
time_it <- function(expr, reps = 5) {
  expr <- substitute(expr)
  env <- parent.frame()
  stats::median(replicate(reps, system.time(eval(expr, env))[["elapsed"]]))
}

compare_nb <- function(d) {
  ours <- nbreg(d$nb, data = d$data)
  ref <- glm.nb(d$nb, data = d$data,
                control = glm.control(epsilon = 1e-10, maxit = 100))
  est <- data.frame(
    term = names(coef(ours)),
    nbreg = coef(ours), glm.nb = coef(ref),
    se_nbreg = sqrt(diag(vcov(ours))), se_glm.nb = sqrt(diag(vcov(ref))),
    row.names = NULL
  )
  est <- rbind(est, data.frame(term = "theta", nbreg = ours$theta,
                               glm.nb = ref$theta, se_nbreg = ours$SE.theta,
                               se_glm.nb = ref$SE.theta))
  diffs <- c(
    coefficients = max(abs(coef(ours) - coef(ref))),
    std_errors = max(abs(sqrt(diag(vcov(ours))) - sqrt(diag(vcov(ref))))),
    theta = abs(ours$theta - ref$theta),
    loglik = abs(as.numeric(logLik(ours)) - as.numeric(logLik(ref)))
  )
  timing <- c(
    ours = time_it(nbreg(d$nb, data = d$data)),
    reference = time_it(glm.nb(d$nb, data = d$data,
                               control = glm.control(epsilon = 1e-10,
                                                     maxit = 100)))
  )
  list(estimates = est, max_abs_diff = diffs, timing = timing,
       loglik = c(ours = as.numeric(logLik(ours)),
                  reference = as.numeric(logLik(ref))),
       iterations = ours$iter)
}

compare_zinb <- function(d) {
  ours <- zinbreg(d$zinb, data = d$data)
  ref <- zeroinfl(d$zinb, data = d$data, dist = "negbin",
                  control = zeroinfl.control(reltol = 1e-12))
  cf <- coef(ours)
  se_ours <- sqrt(diag(vcov(ours)))
  se_ref <- sqrt(diag(vcov(ref)))[names(cf)]
  est <- data.frame(term = names(cf), zinbreg = cf,
                    zeroinfl = coef(ref)[names(cf)],
                    se_zinbreg = se_ours, se_zeroinfl = se_ref,
                    row.names = NULL)
  est <- rbind(est, data.frame(term = "Log(theta)", zinbreg = log(ours$theta),
                               zeroinfl = log(ref$theta),
                               se_zinbreg = ours$SE.logtheta,
                               se_zeroinfl = ref$SE.logtheta))
  diffs <- c(
    coefficients = max(abs(cf - coef(ref)[names(cf)])),
    std_errors = max(abs(c(se_ours - se_ref,
                           ours$SE.logtheta - ref$SE.logtheta))),
    log_theta = abs(log(ours$theta) - log(ref$theta)),
    loglik = abs(as.numeric(logLik(ours)) - as.numeric(logLik(ref)))
  )
  timing <- c(
    ours = time_it(zinbreg(d$zinb, data = d$data), reps = 3),
    reference = time_it(zeroinfl(d$zinb, data = d$data, dist = "negbin",
                                 control = zeroinfl.control(reltol = 1e-12)),
                        reps = 3)
  )
  list(estimates = est, max_abs_diff = diffs, timing = timing,
       loglik = c(ours = as.numeric(logLik(ours)),
                  reference = as.numeric(logLik(ref))),
       iterations = ours$iter)
}

results <- lapply(datasets, function(d) {
  list(nb = compare_nb(d), zinb = compare_zinb(d))
})

for (nm in names(results)) {
  cat("\n==", nm, "==\n")
  cat("NB   max |diff|:\n"); print(signif(results[[nm]]$nb$max_abs_diff, 3))
  cat("ZINB max |diff|:\n"); print(signif(results[[nm]]$zinb$max_abs_diff, 3))
  cat("Timing (s), NB:", results[[nm]]$nb$timing,
      " ZINB:", results[[nm]]$zinb$timing, "\n")
}

saveRDS(results, file.path("inst", "benchmark", "benchmark_results.rds"))
