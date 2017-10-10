library(testthat)
library(ParallelPruning)
library(microbenchmark)


context("PruneTree")
set.seed(1)

EPS <- 1e-8
N <- 2000
tree <- ape::rtree(N)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rexp(N, 1/.01)
z <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sqrt(sigmae^2+se^2))

pruneInfo2 <- POUMM::pruneTree(tree, z, se)
pruneInfo <- ParallelPruning::prepare(tree)

ppa <- ParallelPruning:::ParallelPruningAlgorithm$new()

poummLikelihood <- ParallelPruning:::POUMM_Likelihood$new()

poummLikelihood$set_treeAndData(tree = tree, z = z[1:N], se = se[1:N])

poummLikelihood$do_pruning(
  list(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae))


expect_lt(abs(
  POUMM:::loglik_abc_g0_g0Prior(
    abc = poummLikelihood$abc(), g0Prior = NA,
    g0 = g0, alpha = alpha, theta = theta, sigma = sigma)$loglik -
    POUMM::likPOUMMGivenTreeVTips(
      z = z, tree = tree,
      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sqrt(sigmae^2+se^2),
      g0Prior = NA
    )), EPS)


#
# all(pruneInfo2$integrator$eReord-pruneInfo$eReord == -1)
#
# all(pruneInfo2$implementationCPP$tReord==pruneInfo$tReord)
#
# all(pruneInfo2$implementationCPP$zReord[1:N]==pruneInfo$z[1:N][pruneInfo2$reord])
#
# all(pruneInfo2$implementationCPP$seReord==pruneInfo$seReord[pruneInfo2$reord])
#
