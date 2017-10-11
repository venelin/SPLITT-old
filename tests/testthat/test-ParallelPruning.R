library(testthat)
library(ParallelPruning)
library(microbenchmark)


lik_POUMM <- function(
  poummObj,
  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae) {

  poummObj$do_pruning(
    list(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae))

  POUMM:::loglik_abc_g0_g0Prior(
    abc = poummObj$abc(), g0Prior = NA,
    g0 = g0, alpha = alpha, theta = theta, sigma = sigma)$loglik
}

lik_POUMM_old <- function(
  pruneInfo,
  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae) {

  POUMM::likPOUMMGivenTreeVTipsC4(
    integrator = pruneInfo$integrator,
    g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae)
}

context("PruneTree")
set.seed(1)

EPS <- 1e-8
N <- 10
tree <- ape::rtree(N)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rexp(N, 1/.01)
z <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)


microbenchmark(
  pruneInfo <- POUMM::pruneTree(tree, z, se),
  poummLikelihood <- ParallelPruning:::POUMM_Likelihood$new(tree, z = z[1:N], se = se),
  poummabc <- ParallelPruning:::POUMM_abc$new(tree, z = z[1:N], se = se),
  times = 1)

# test correct value
test_that("POUMM likelihood", expect_lt(abs(
  lik_POUMM(poummLikelihood,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))

test_that("POUMM abc", expect_lt(abs(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))

microbenchmark(
  lik_POUMM(poummLikelihood,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),
  lik_POUMM_old(pruneInfo,
                g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae))
