bio_formula <- art ~ fem + mar + kid5 + phd + ment

test_that("nbreg matches MASS::glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- nbreg(bio_formula, data = bioChemists)
  ref <- MASS::glm.nb(bio_formula, data = bioChemists,
                      control = stats::glm.control(epsilon = 1e-10, maxit = 100))
  expect_true(fit$converged)
  expect_lte(max(abs(coef(fit) - coef(ref))), 1e-5)
  expect_lte(max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))), 1e-5)
  expect_lte(abs(fit$theta - ref$theta), 1e-4)
  expect_lte(abs(fit$SE.theta - ref$SE.theta), 1e-4)
  expect_lte(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
  expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
  expect_equal(AIC(fit), AIC(ref), tolerance = 1e-6)
  expect_equal(BIC(fit), BIC(ref), tolerance = 1e-6)
})

test_that("nbreg matches MASS::glm.nb on NMES1988", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER")
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  fit <- nbreg(f, data = NMES1988)
  ref <- MASS::glm.nb(f, data = NMES1988,
                      control = stats::glm.control(epsilon = 1e-10, maxit = 100))
  expect_lte(max(abs(coef(fit) - coef(ref))), 1e-5)
  expect_lte(max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))), 1e-5)
  expect_lte(abs(fit$theta - ref$theta), 1e-4)
  expect_lte(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
})

test_that("nbreg handles offsets in the formula and as an argument", {
  skip_if_not_installed("MASS")
  set.seed(651)
  n <- 400
  dat <- data.frame(x = rnorm(n), expo = runif(n, 0.5, 3))
  dat$y <- rnbinom(n, size = 2, mu = dat$expo * exp(0.3 + 0.6 * dat$x))
  f1 <- nbreg(y ~ x + offset(log(expo)), data = dat)
  f2 <- nbreg(y ~ x, data = dat, offset = log(expo))
  ref <- MASS::glm.nb(y ~ x + offset(log(expo)), data = dat,
                      control = stats::glm.control(epsilon = 1e-10))
  expect_equal(coef(f1), coef(f2), tolerance = 1e-8)
  expect_lte(max(abs(coef(f1) - coef(ref))), 1e-5)
  # predictions use the offset from newdata
  nd <- data.frame(x = c(0, 1), expo = c(1, 2))
  expect_equal(unname(predict(f1, nd, type = "response")),
               unname(predict(ref, nd, type = "response")), tolerance = 1e-5)
  expect_equal(predict(f2, nd, type = "response"),
               predict(f1, nd, type = "response"), tolerance = 1e-8)
})

test_that("nbreg works with a single covariate and with prior weights", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- nbreg(art ~ ment, data = bioChemists)
  ref <- MASS::glm.nb(art ~ ment, data = bioChemists,
                      control = stats::glm.control(epsilon = 1e-10))
  expect_lte(max(abs(coef(fit) - coef(ref))), 1e-5)

  w <- rep(1:3, length.out = nrow(bioChemists))
  fw <- nbreg(art ~ fem + ment, data = bioChemists, weights = w)
  rw <- MASS::glm.nb(art ~ fem + ment, data = bioChemists, weights = w,
                     control = stats::glm.control(epsilon = 1e-10))
  expect_lte(max(abs(coef(fw) - coef(rw))), 1e-5)
  expect_lte(abs(fw$theta - rw$theta), 1e-4)
})

test_that("nbreg gives informative errors for bad input", {
  dat <- data.frame(x = 1:10, y = c(0:8, -1))
  expect_error(nbreg(y ~ x, data = dat), "non-negative")
  dat$y <- c(0:8, 2.5)
  expect_error(nbreg(y ~ x, data = dat), "integer")
  dat$y <- 0:9
  expect_error(nbreg(y ~ x, data = dat, weights = c(-1, rep(1, 9))),
               "negative weights")
  dat$x2 <- 2 * dat$x
  expect_error(nbreg(y ~ x + x2, data = dat), "rank deficient")
  expect_error(nbreg_control(epsilon = -1), "epsilon")
})

test_that("nbreg warns in the Poisson limit and when maxit is reached", {
  set.seed(1)
  x <- rnorm(300)
  y <- rpois(300, exp(1 + 0.3 * x))
  expect_warning(nbreg(y ~ x, control = nbreg_control(theta_max = 5)),
                 "theta")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  expect_warning(nbreg(art ~ ment, data = bioChemists,
                       control = nbreg_control(maxit = 1)),
                 "did not converge")
})

test_that("nbreg methods return consistent results", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- nbreg(bio_formula, data = bioChemists)
  ref <- MASS::glm.nb(bio_formula, data = bioChemists)
  expect_equal(nobs(fit), nrow(bioChemists))
  expect_equal(unname(fitted(fit)), unname(fitted(ref)), tolerance = 1e-5)
  expect_equal(unname(residuals(fit, "pearson")),
               unname(residuals(ref, "pearson")), tolerance = 1e-5)
  expect_equal(unname(residuals(fit, "deviance")),
               unname(residuals(ref, "deviance")), tolerance = 1e-5)
  expect_equal(residuals(fit, "response"), bioChemists$art - fitted(fit))
  expect_equal(predict(fit, type = "response"), fitted(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_output(print(fit), "Theta")
  expect_output(print(s), "Coefficients")
})
