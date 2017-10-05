library(testthat)
library(ParallelPruning)
library(microbenchmark)


context("PruneTree")
set.seed(1)

N <- 10000
tree <- ape::rtree(N)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rexp(N, 1/.01)
z <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sqrt(sigmae^2+se^2))


microbenchmark(pruneInfo2 <- ParallelPruning::pruneTree2(tree, z, se),
pruneInfo <- ParallelPruning::pruneTree(tree, z, se), times=1)

all(pruneInfo2$implementationCPP$eReord-pruneInfo$eReord == -1)

all(pruneInfo2$implementationCPP$tReord==pruneInfo$tReord)

all(pruneInfo2$implementationCPP$zReord[1:N]==pruneInfo$z[1:N][pruneInfo2$reord])

all(pruneInfo2$implementationCPP$seReord==pruneInfo$seReord[pruneInfo2$reord])

