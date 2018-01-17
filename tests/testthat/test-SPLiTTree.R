library(testthat)
library(SPLiTTree)
library(microbenchmark)
library(ape)
library(POUMM)

lik_POUMM <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  abc = poummObj$TraverseTree(c(alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae), mode)

  POUMM:::loglik_abc_g0_g0Prior(
    abc = abc, g0Prior = NA,
    g0 = g0, alpha = alpha, theta = theta, sigma = sigma)$loglik
}

lik_POUMM_lnDetV_Q <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  res <- poummObj$TraverseTree(c(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),mode)

  -1/2*(poummObj$tree$num_tips*log(2*pi) + 2*res[3] + res[1]+res[2])
}

lik_POUMM_old <- function(
  pruneInfo,
  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae) {

  POUMM::likPOUMMGivenTreeVTipsC(
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
z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)

# load("../../../poummBayesianValidation/DATA/Microbenchmark/Trees.RData")
# tree<-trees$tree[[24]]
# z<-trees$z[[24]]


#ppt<-SPLiTTree:::ParallelPruningTree$new(tree)
context("POUMM::pruneTree")
pruneInfo <- POUMM::pruneTree(tree, z, se)


context("SPLiTTree:::ParallelPruningThreePointPOUMM")
poummlnDetVQ <- SPLiTTree:::ParallelPruningThreePointPOUMM$new(tree, z = z[1:N], se = se)

# test correct value
test_that("POUMM lnDetV_Q", expect_lt(abs(
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
            g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21) -
    lik_POUMM_old(pruneInfo,
                  g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
    )), EPS))


test_that("Equal parallel vs serial pruning", expect_equal(
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
                     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21),
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
                     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 1)))


library(ape)
library(SPLiTTree)
library(microbenchmark)
library(testthat)
N <- 10000
tree<-rtree(N)
idx <- sample(which(tree$edge[, 2] > N), size=N/2)
tree$edge.length[idx] <- 0
#tree2 <- di2multi(tree)
z <- rnorm(N)
se <- rep(0, N)
context("PruneTree")
set.seed(1)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7

lik_POUMM_lnDetV_Q <- function(
  poummObj,
  g0, alpha, theta, sigma, sigmae, mode) {

  res <- poummObj$TraverseTree(c(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),mode)

  -1/2*(N*log(2*pi) + 2*res[3] + res[1]+res[2])
}

poummlnDetVQ <- SPLiTTree:::ParallelPruningThreePointPOUMM$new(tree, z = z[1:N], se = se)



microbenchmark(
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
                     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 11),
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
                     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 12),
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
                     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21),
  lik_POUMM_lnDetV_Q(poummlnDetVQ,
                     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 22))
