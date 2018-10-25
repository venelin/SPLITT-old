## ----setup, include = FALSE-------------------------------------------------------------------------------------------
# Make results reproducible
set.seed(1)
knitr::opts_chunk$set(cache = FALSE)
options(digits = 4, width=120, scipen = 999)
remove.packages("SPLITT")

## ---- echo=TRUE, eval=FALSE-------------------------------------------------------------------------------------------
#  devtools::install_github("venelin/SPLITT")

## ---- echo=FALSE, eval=TRUE-------------------------------------------------------------------------------------------
capture <- capture.output( devtools::install_github("venelin/SPLITT") )
capture <- 'Building SPLITT vignettes
* installing *source* package ‘SPLITT’ ...
** libs
/usr/local/clang6/bin/clang++  -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -fopenmp -I"/Users/vmitov/Library/R/3.5/library/Rcpp/include" -I"/Users/vmitov/Library/R/3.5/library/RcppArmadillo/include" -I/usr/local/include   -fPIC  -Wall -g -O2  -c RCPP__AbcPMM.cpp -o RCPP__AbcPMM.o
/usr/local/clang6/bin/clang++  -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -fopenmp -I"/Users/vmitov/Library/R/3.5/library/Rcpp/include" -I"/Users/vmitov/Library/R/3.5/library/RcppArmadillo/include" -I/usr/local/include   -fPIC  -Wall -g -O2  -c RCPP__Tree.cpp -o RCPP__Tree.o
/usr/local/clang6/bin/clang++  -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -fopenmp -I"/Users/vmitov/Library/R/3.5/library/Rcpp/include" -I"/Users/vmitov/Library/R/3.5/library/RcppArmadillo/include" -I/usr/local/include   -fPIC  -Wall -g -O2  -c RcppExports.cpp -o RcppExports.o
/usr/local/clang6/bin/clang++ -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o SPLITT.so RCPP__AbcPMM.o RCPP__Tree.o RcppExports.o -fopenmp -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
installing to /Users/vmitov/Library/R/3.5/library/SPLITT/libs
** R
** tests
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded
* DONE (SPLITT)'

cat(capture)

## ---------------------------------------------------------------------------------------------------------------------
library(SPLITT)
MiniBenchmark(N=10, Ntests = 1000)

## ----treeAndDataLogLik------------------------------------------------------------------------------------------------
library(ape)
newick <- '(((1:3,2:2.8)6:4.2,3:4.1)7:4.5,(4:4,5:5)8:6)0;'
tree <- read.tree(text=newick)

x0 <- 0.1
sigma2 <- 0.25
sigmae2 <- 1

set.seed(1)
g <- rTraitCont(tree, model = "BM", root.value = x0, 
                sigma = sqrt(sigma2),
                ancestor = FALSE)

x <- g + rnorm(n = length(tree$tip.label), mean = 0, sd = sqrt(sigmae2))

PMMLogLikCpp(x, tree, x0, sigma2, sigmae2, mode = 1)

## ----fig1a, dpi=150, fig.width=4.6------------------------------------------------------------------------------------
par(mfrow=c(1,5))
par(mar=c(0,0,0,0))
par(oma=c(0, 0,0, 0))

PlotParallelTraversal(tree)

