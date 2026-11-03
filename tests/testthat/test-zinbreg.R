bio_zformula <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment

compare_zinb <- function(fit, ref) {
  cf <- coef(fit)
  expect_lte(max(abs(cf - coef(ref)[names(cf)])), 1e-3)
  expect_lte(max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref)))[names(cf)])),
             1e-3)
  expect_lte(abs(log(fit$theta) - log(ref$theta)), 1e-3)
  expect_lte(abs(fit$SE.logtheta - ref$SE.logtheta), 1e-3)
  expect_lte(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-3)
}

test_that("zinbreg matches pscl::zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(bio_zformula, data = bioChemists)
  ref <- pscl::zeroinfl(bio_zformula, data = bioChemists, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  expect_true(fit$converged)
  compare_zinb(fit, ref)
  expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
})

test_that("zinbreg matches pscl::zeroinfl on NMES1988", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER")
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  fit <- zinbreg(f, data = NMES1988)
  ref <- pscl::zeroinfl(f, data = NMES1988, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  compare_zinb(fit, ref)
})

test_that("EM never decreases the observed-data log-likelihood", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(bio_zformula, data = bioChemists)
  expect_true(fit$monotone)
  expect_true(all(diff(fit$loglik_trace) >= -1e-8))
})

test_that("analytic Hessian agrees with a numerical Hessian", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + kid5 + ment | ment + kid5, data = bioChemists)
  par <- c(fit$coefficients$count, fit$coefficients$zero, log(fit$theta))
  Ha <- zinb_hessian(fit$x, fit$z, fit$y, fit$offset, fit$weights,
                     fit$coefficients$count, fit$coefficients$zero, fit$theta)
  Hn <- numDeriv::hessian(zinb_loglik_par, par, X = fit$x, Z = fit$z,
                          y = fit$y, offset = fit$offset, a = fit$weights)
  expect_lte(max(abs(Ha - Hn)) / max(abs(Hn)), 1e-5)
  fit_num <- zinbreg(art ~ fem + kid5 + ment | ment + kid5,
                     data = bioChemists,
                     control = zinbreg_control(hessian = "numerical"))
  expect_equal(sqrt(diag(vcov(fit_num))), sqrt(diag(vcov(fit))),
               tolerance = 1e-4)
})

test_that("a formula without '|' uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  f1 <- zinbreg(art ~ kid5 + ment, data = bioChemists)
  f2 <- zinbreg(art ~ kid5 + ment | kid5 + ment, data = bioChemists)
  expect_equal(names(coef(f1, "zero")), c("(Intercept)", "kid5", "ment"))
  expect_equal(coef(f1), coef(f2))
})

test_that("zinbreg handles a count-part offset and single covariates", {
  skip_if_not_installed("pscl")
  set.seed(651)
  n <- 500
  dat <- data.frame(x = rnorm(n), w = rnorm(n), expo = runif(n, 0.5, 2))
  pi <- plogis(-1 + 0.8 * dat$w)
  mu <- dat$expo * exp(1 + 0.5 * dat$x)
  dat$y <- ifelse(runif(n) < pi, 0, rnbinom(n, size = 1.5, mu = mu))
  fit <- zinbreg(y ~ x + offset(log(expo)) | w, data = dat)
  ref <- pscl::zeroinfl(y ~ x + offset(log(expo)) | w, data = dat,
                        dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  compare_zinb(fit, ref)
  fit2 <- zinbreg(y ~ x | w, data = dat, offset = log(expo))
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)
  nd <- data.frame(x = c(0, 1), w = c(0, 1), expo = c(1, 2))
  for (type in c("response", "count", "zero")) {
    expect_equal(unname(predict(fit, nd, type = type)),
                 unname(predict(ref, nd, type = type)), tolerance = 1e-3)
  }
})

test_that("zinbreg gives informative errors for bad input", {
  dat <- data.frame(x = rnorm(20), y = rpois(20, 3) + 1)
  expect_error(zinbreg(y ~ x, data = dat), "no zeros")
  dat$y[1] <- -1
  expect_error(zinbreg(y ~ x, data = dat), "non-negative")
  dat$y[1] <- 0
  expect_error(zinbreg(y ~ x | x | x, data = dat), "at most two parts")
  expect_error(zinbreg_control(maxit = 0), "maxit")
})

test_that("zinbreg warns when EM hits maxit", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  expect_warning(zinbreg(art ~ ment, data = bioChemists,
                         control = zinbreg_control(maxit = 2)),
                 "did not converge")
})

test_that("zinbreg methods agree with zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(bio_zformula, data = bioChemists)
  ref <- pscl::zeroinfl(bio_zformula, data = bioChemists, dist = "negbin")
  expect_equal(nobs(fit), nrow(bioChemists))
  expect_equal(unname(fitted(fit)), unname(fitted(ref)), tolerance = 1e-4)
  expect_equal(unname(residuals(fit)), unname(residuals(ref)), tolerance = 1e-4)
  expect_equal(unname(residuals(fit, "response")),
               unname(residuals(ref, "response")), tolerance = 1e-4)
  expect_equal(predict(fit, type = "response"), fitted(fit))
  expect_equal(unname(predict(fit, type = "zero")),
               unname(predict(ref, type = "zero")), tolerance = 1e-4)
  expect_equal(unname(predict(fit, type = "prob")),
               unname(predict(ref, type = "prob")), tolerance = 1e-4)
  expect_equal(dim(vcov(fit, "count")), c(6L, 6L))
  expect_equal(names(coef(fit, "zero")), rownames(vcov(fit, "zero")))
  s <- summary(fit)
  expect_s3_class(s, "summary.zinbreg")
  expect_equal(rownames(s$coefficients$count)[7], "Log(theta)")
  expect_output(print(fit), "Zero-inflation")
  expect_output(print(s), "Log\\(theta\\)")
})
