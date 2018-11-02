#library(RhpcBLASctl)
library(data.table)
#library(Rphylopars)
library(data.table)
# library(PCMBase)
# library(PCMBaseCpp)
source("timeXXX.R")

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  args <- as.integer(args)
  cat("args = ", toString(args))
} else {
  # args[1] is the number of openmp threads; args[2] is the number of blas threads.
  args <- c(omp_get_max_threads(), 1)
}

saveResults <- TRUE

if(R.version[['os']]=='linux-gnu') {
  # this only works on linux
  cpuInfo <- system("cat /proc/cpuinfo | grep 'model name' | uniq", intern = TRUE)
} else {
  # this only works on mac OS x
  cpuInfo <- system("sysctl -a -n machdep.cpu.brand_string", intern = TRUE)
}

cat("Process started with settings: \n")
cat(cpuInfo, "\n")
cat("num_cores = ", get_num_cores(), "\n")
cat("num_procs = ", get_num_procs(), "\n")
cat("blas_num_procs = ", blas_get_num_procs(), "\n")
cat("omp_num_procs = ", omp_get_num_procs(), "\n")
cat("omp_max_threads = ", omp_get_max_threads(), "\n")

if(args[1] != omp_get_max_threads()) {
  cat("omp_set_num_threads to:", args[1], "\n")
  omp_set_num_threads(args[1])
}
if(args[2] != 0) {
  cat("blas_set_num_threads to:", args[2], "\n")
  blas_set_num_threads(args[2])
}

nTests <- 5 # number of executions in each microbenchmark

times <- values <- trees <- NULL

# load the trees data.table
load('Trees.RData')
if(length(args) >= 3) {
  ids <- args[-(1:2)]
} else {
  ids <- 1:nrow(trees)
}


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
#x <- timePCMBaseCpp(X, tree, model, 11)
#x <- timePCMBaseCpp(X, tree, model, 21)
rm(x)

gc()

for(id in ids) {
  tree <- trees$tree[[id]]

  X <- trees$X[[id]]
  N <- trees$N[[id]]
  k <- trees$num_traits[[id]]
  #model <- trees$model[[id]]

  #print(model)
  collessNorm <- trees$collessNorm[[id]]

  pSymb <- trees$pSymb[[id]]

  cat("Tree-id:", id, "; tree-size:", trees$N[[id]], " Colless:", trees$collessNorm[[id]], ", num_traits=", k, "\n")

  if(k == 1 && 
     (id != 13 || N < 1e5))
  times <- rbind(
    times,
    cbind(data.table(id = id,
                     cpuInfo=cpuInfo,
                     num_cores = get_num_cores(),
                     num_procs = get_num_procs(),
                     blas_num_procs = blas_get_num_procs(),
                     omp_num_procs = omp_get_num_procs(),
                     omp_max_threads = omp_get_max_threads(),
                     pSymb = pSymb, collessNorm = collessNorm, N = N, k = k),
          rbindlist(list(
            if(k == 1) timePOUMMCpp(X, tree, model, 11) else NULL,
            if(k == 1) timePOUMMCpp(X, tree, model, 21) else NULL,
            if(k == 1) timePOUMMCpp(X, tree, model, 24) else NULL,
            if(k == 1) timePOUMMCpp(X, tree, model, 25) else NULL,
            if(k == 1) timePOUMMCpp(X, tree, model, 31) else NULL

            # timePCMBaseCpp(X, tree, model, 11),
            # timePCMBaseCpp(X, tree, model, 21),
            # timePCMBaseCpp(X, tree, model, 24),
            # timePCMBaseCpp(X, tree, model, 25),
            # timePCMBaseCpp(X, tree, model, 31)
            # timeRphylopars(X, tree)
          ))))

  print(times)

  if(saveResults) {
    save(times,
         file = paste0("Result-12-cluster-omp-", args[1], "-blas-", args[2], ".RData"))
  }
}
