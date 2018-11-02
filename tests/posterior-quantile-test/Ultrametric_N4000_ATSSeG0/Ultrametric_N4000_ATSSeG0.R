# source this file within the simName directory;
# then execute code in the if(FALSE) chunks manually.

library(POUMM)
library(data.table)

# this script must be executed under directory DATA/*simName*
simName <- "Ultrametric_N4000_ATSSeG0"
treeFileName <- "tree_bd_taxa_N4000_2_1_1.RData"

if(FALSE) {
  setwd(paste0("~/Documents/Bio/Projects/poummBayesianValidation/DATA/", simName))
}

if(FALSE) {
  # this code block is executed manually in order to generate the tree
  library(TreeSim)
  set.seed(1)
  
  # 1.
  # Number of tips
  N <- 4000
  tree <- TreeSim::sim.bd.taxa(n = N, numbsim = 1, lambda = 2, mu = 1, frac = 1,
                               complete = FALSE)[[1]]
  
  treeFilePath <- paste0(treeFileName)
  save(tree, file = treeFilePath)
}

load(treeFileName)
tMean <- mean(nodeTimes(tree, tipsOnly = TRUE))

if(FALSE) {
  # This is a manually executed code-block in order to generate the benchmark
  # and transfer the files to the hpc.
  nReplics <- 2000
  
  data <- data.table(id = 1:nReplics, key="id")
  
  save(data, file = paste0(simName, ".RData"))
  
  ids <- data[, id]
  
  genBenchJobs(replication, script.file = simName, 
               table.file = paste0(simName, ".RData"), 
               ids = ids, 
               perJob = 1, 
               type = "bsub", 
               requires = c("POUMM"), 
               sources = c(paste0(simName, ".R"), 
                           "../../defineSimPMMReplication_ATSSeG0.R"),
               bsub.W = "24:00",
               bsub.mem = 6000,
               bsub.other = "-R beta ",
               sleep.every = 200)
}

execName <- "exec_PMM_20170428"

if(FALSE) {
  # transfer files to euler
  #system('rsync -cv -r --progress ~/Documents/Bio/Projects/poumm_1.0.0.tar.gz vmitov@euler:~/')
  #system('rsync -cv -r --progress ~/Documents/Bio/Projects/patherit_0.7.tar.gz vmitov@euler:~/')
  #system('rsync -cv -r --progress ~/Documents/Bio/Projects/benchtable_1.0.tar.gz vmitov@euler:~/')
  
  system(paste0('rsync -cv -r --progress ', '../defineSimReplication_ATSSeG0.R',  ' vmitov@euler:~/poummBayesianValidation/DATA2/'))
  
  system(paste0('rsync -cv -r --progress ', simName, '.R',  ' vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  system(paste0('rsync -cv -r --progress ', simName,  '.RData vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  system(paste0('rsync -cv -r --progress ', treeFileName,  ' vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  system(paste0('rsync -cv -r --progress j_* vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  
} 

if(FALSE) {
  # bring only job files from euler 
  system(paste0('rsync -cv -r --remove-source-files --progress vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/job* ', './', execName))
  
  # bring all files in simulation execution
  system(paste0('rsync -cv -r --remove-source-files --progress vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/* ', './', execName))
}


