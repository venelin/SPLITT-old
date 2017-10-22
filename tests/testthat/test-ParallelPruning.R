library(testthat)
library(ParallelPruning)
library(microbenchmark)


lik_POUMM <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  poummObj$set_parameters(alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae)
  poummObj$do_pruning(mode)

  POUMM:::loglik_abc_g0_g0Prior(
    abc = poummObj$abc(), g0Prior = NA,
    g0 = g0, alpha = alpha, theta = theta, sigma = sigma)$loglik
}

lik_POUMM_lnDetV_Q <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  poummObj$set_parameters(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae)
  poummObj$do_pruning(mode)

  -1/2*(poummObj$N*log(2*pi) + 2*poummObj$lnDetD + poummObj$lnDetV+poummObj$Q)
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
N <- 10000
tree <- ape::rtree(N)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rep(0, N); #rexp(N, 1/.01)
z <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)
pruneInfo <- POUMM::pruneTree(tree, z, se)

poummabc <- ParallelPruning:::POUMM_abc$new(tree, z = z[1:N], se = se)

poummabc$set_parameters( alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae)

poummabc$do_pruning(0)
poummabc$abc()

poummlnDetVQ <- ParallelPruning:::POUMM_lnDetV_Q_1d$new(tree, z = z[1:N])

# test correct value
test_that("POUMM abc", expect_lt(abs(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 2) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))

test_that("POUMM lnDetV_Q", expect_lt(abs(
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 2) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))


test_that("Equal parallel vs serial pruning", expect_equal(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode =2),
    lik_POUMM(poummabc,
              g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 1)))

test_that("Equal hybrid vs serial pruning", expect_equal(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode =0),
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 1)))



microbenchmark(
  pruneInfo <- POUMM::pruneTree(tree, z, se),
  times = 1)


microbenchmark(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 1),
  lik_POUMM_old(pruneInfo,
                g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae))
