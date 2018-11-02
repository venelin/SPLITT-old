#library(diversitree)
#library(geiger)
library(data.table)
# library(phylolm)
# library(Rphylopars)
# library(PCMBase)
# library(PCMBaseCpp)
source("timeXXX.R")

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  ids <- as.integer(args)
}

compilerInfo <- "icpc-omp-for-simd"
saveResults <- TRUE

if(R.version[['os']]=='linux-gnu') {
  # this only works on linux
  cpuInfo <- system("cat /proc/cpuinfo | grep 'model name' | uniq", intern = TRUE)
} else {
  # this only works on mac OS x
  cpuInfo <- system("sysctl -a -n machdep.cpu.brand_string", intern = TRUE)
}

print(cpuInfo)

nTests <- 5 # number of executions in each microbenchmark

times <- values <- trees <- NULL

# load the trees data.table
load('Trees.RData')

# needed for warmup -loading of C++ libraries
cat("Warming up on id = 1 ")
id <- 1

tree <- trees$tree[[id]]

X <- trees$X[[id]]
N <- trees$N[[id]]
k <- trees$num_traits[[id]]
model <- trees$model[[id]]

print(model)

collessNorm <- trees$collessNorm[[id]]

pSymb <- trees$pSymb[[id]]

x <- timePOUMMR(X, tree, model)
#x <- timePCMBaseR(X, tree, model)
x <- timePOUMMCpp(X, tree, model, 11)
x <- timePOUMMCpp(X, tree, model, 21)
#x <- timeRphylopars(X=X, tree = tree)
#x <- timePCMBaseCpp(X, tree, model, 11)
#x <- timePCMBaseCpp(X, tree, model, 21)
rm(x)

gc()

for(id in ids) {
  tree <- trees$tree[[id]]

  X <- trees$X[[id]]
  N <- trees$N[[id]]
  k <- trees$num_traits[[id]]
  model <- trees$model[[id]]

  print(model)
  collessNorm <- trees$collessNorm[[id]]

  pSymb <- trees$pSymb[[id]]

  cat("Tree-id:", id, "; tree-size:", trees$N[[id]], " Colless:", trees$collessNorm[[id]], ", num_traits=", k, "\n")

  times <- rbind(
    times,
    cbind(data.table(compilerInfo = compilerInfo,
                     cpuInfo=cpuInfo, nCores = 4, pSymb = pSymb, collessNorm = collessNorm, N = N, k = k),
          rbindlist(list(
            timePOUMMR(X, tree, model),
            #timePCMBaseR(X, tree, model),
            timePOUMMCpp(X, tree, model, 11),
            timePOUMMCpp(X, tree, model, 21),
            timePOUMMCpp(X, tree, model, 25)
            # timeRphylopars(X=X, tree = tree),
            # timePCMBaseCpp(X, tree, model, 11),
            # timePCMBaseCpp(X, tree, model, 21),
            # timePCMBaseCpp(X, tree, model, 25)
          ))))

  print(times)
}

if(saveResults) {
  times1 <- copy(times)

  if(!any(ids == 1)) {
    load("Results-microbenchmark-Mac-1.RData")
  }

  times <- rbind(times, times1)
  save(times,
       file = paste0("Results-microbenchmark-Mac-1.RData"))
}
