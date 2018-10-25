library(ape)
library(SPLITT)
library(POUMM)

set.seed(10)

N <- 100
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

POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2)
POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2)

if(require(POUMM)) {
  likPOUMMGivenTreeVTips(x, tree, alpha, theta, sqrt(sigma2), sqrt(sigmae2), g0 = x0)
}
  