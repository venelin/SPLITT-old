library(testthat)
library(ParallelPruning)
library(microbenchmark)
library(ape)

lik_POUMM <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  abc = poummObj$DoPruning(c(alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae), mode)

  POUMM:::loglik_abc_g0_g0Prior(
    abc = abc, g0Prior = NA,
    g0 = g0, alpha = alpha, theta = theta, sigma = sigma)$loglik
}

lik_POUMM_lnDetV_Q <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  res <- poummObj$DoPruning(c(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),mode)

  -1/2*(poummObj$tree$num_tips*log(2*pi) + 2*res[3] + res[1]+res[2])
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
N <- 100
tree <- ape::rtree(N)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rep(0, N); #rexp(N, 1/.01)
z <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)

# load("../../../poummBayesianValidation/DATA/Microbenchmark/Trees.RData")
# tree<-trees$tree[[24]]
# z<-trees$z[[24]]


#ppt<-ParallelPruning:::ParallelPruningTree$new(tree)

pruneInfo <- POUMM::pruneTree(tree, z, se)

poummabc <- ParallelPruning:::ParallelPruningAbcPOUMM$new(tree, z = z[1:N], se = se)

poummlnDetVQ <- ParallelPruning:::ParallelPruningThreePointPOUMM$new(tree, z = z[1:N], se = se)

# test correct value
test_that("POUMM abc", expect_lt(abs(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))

test_that("POUMM lnDetV_Q", expect_lt(abs(
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))


test_that("Equal parallel vs serial pruning", expect_equal(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode =21),
    lik_POUMM(poummabc,
              g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 10)))

test_that("Equal hybrid vs serial pruning", expect_equal(
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode =31),
  lik_POUMM(poummabc,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 12)))
