context("wevid")

data(cleveland)
data(pima)
data(fitonly)

source("dump_wevid.R")

test_that("linear regression",
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
