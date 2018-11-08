library(testthat)
context("Test R and Cpp code calculate the same PMM log-likelihood")

library(ape)
library(PMMUsingSPLITT)

set.seed(10)

N <- 1000
x0 <- 0.1
alpha <- 1
theta <- 10
sigma2 <- 0.25
sigmae2 <- 1

tree <- rtree(N)

g <- rTraitCont(tree, model = "OU", root.value = x0,
                alpha = alpha, sigma = sqrt(sigma2),
                ancestor = FALSE)

x <- g + rnorm(n = N, mean = 0, sd = sqrt(sigmae2))

test_that(
  "PMMLogLik == PMMLogLikCpp",
  expect_equal(PMMLogLik(x, tree, x0, sigma2, sigmae2),
               PMMLogLikCpp(x, tree, x0, sigma2, sigmae2))
)
