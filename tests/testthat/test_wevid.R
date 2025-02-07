context("wevid")

data(cleveland)
data(pima)
data(fitonly)

source("dump_wevid.R")

test_that("example datasets",
{
  d1 <- with(cleveland, Wdensities(y, posterior.p, prior.p))
  d2 <- with(pima, Wdensities(y, posterior.p, prior.p))
  d3 <- with(fitonly, Wdensities(y, posterior.p, prior.p))

  s1 <- summary(d1)
  s2 <- summary(d2)
  s3 <- summary(d3)

  expect_equal(s1, s1.ok)
  expect_equal(s2, s2.ok)
  expect_equal(s3, s3.ok)
})

test_that("two components in controls, one in cases",
{
  set.seed(1)
  ex1.y <- c(rep(0, 50), rep(1, 50))
  ex1.posterior.p <- c(runif(50, 0, 0.1), runif(25, 0, 0.2), runif(25, 0.4, 0.6))
  ex1.prior.p <- 0.1
  expect_message(d1 <- Wdensities(ex1.y, ex1.posterior.p, ex1.prior.p),
                 "Density with 2 mixture components chosen by BIC")
  s1 <- summary(d1)

  expect_equal(s1$`Crude C-statistic`, 0.870)
  expect_equal(s1$`Crude Λ (bits)`, 1.38)
  expect_equal(s1$`Test log-likelihood (nats)`, -83.97)
})

test_that("two components in cases and controls, but a component with one point",
{
  set.seed(1)
  y <- c(rep(0, 45), rep(1, 50))
  posterior.p <- c(runif(45, 0, 0.25), runif(50, 0, 0.1))
  prior.p <- 0.1
  expect_message(d1 <- Wdensities(y, posterior.p, prior.p),
                 "Density with 2 mixture components chosen by BIC")
  s1 <- summary(d1)

  expect_equal(s1$`Crude C-statistic`, 0.144)
  expect_equal(s1$`Crude Λ (bits)`, -0.79)
  expect_equal(s1$`Test log-likelihood (nats)`, -168.49)
})

test_that("posterior probabilities 0 or 1",
{
  p10 <- pima
  p10$posterior.p[100] <- 0
  expect_error(with(p10, Wdensities(y, posterior.p, prior.p)),
               "exactly 0 or 1")
  expect_error(with(p10, weightsofevidence(posterior.p, prior.p)),
               "exactly 0 or 1")

  p10$posterior.p[100] <- 1
  expect_error(with(p10, Wdensities(y, posterior.p, prior.p)),
               "exactly 0 or 1")
  expect_error(with(p10, weightsofevidence(posterior.p, prior.p)),
               "exactly 0 or 1")
})
