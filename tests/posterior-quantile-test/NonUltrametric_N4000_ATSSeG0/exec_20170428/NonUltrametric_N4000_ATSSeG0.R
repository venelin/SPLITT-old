# IMPORTANT: when running locally to prepare scripts: 
# 1. frist source this file within the simName directory;
# 2. then execute code in the if(FALSE) chunks manually.

library(POUMM)
library(data.table)

# this script must be executed under directory DATA/*simName*
simName <- "NonUltrametric_N4000_ATSSeG0"
treeFileName <- "tree_bdsky_stt_N4000_2_1_1.RData"

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
  tree <- TreeSim::sim.bdsky.stt(N, lambdasky = 2, deathsky = 1, 
                                 timesky=c(0), sampprobsky = 1)[[1]]
  
  treeFilePath <- paste0(treeFileName)
  save(tree, file = treeFilePath)
}

# excecute this only in jobs where the treeFileName is in the execution directory execName
load(treeFileName)
N <- length(tree$tip.label)
tMean <- mean(nodeTimes(tree, tipsOnly = TRUE))

if(FALSE) {
  # This is a manually executed code-block in order to generate the benchmark
  # and transfer the files to the hpc.
  require(benchtable)
  nReplics <- 2000
  
  data <- data.table(id = 1:nReplics, key="id")
  ids <- data[, id]
  
  save(data, file = paste0(simName, ".RData"))
  
  genBenchJobs("replication", script.file = simName, 
               table.file = paste0(simName, ".RData"), 
               ids = ids, 
               perJob = 1, 
               type = "bsub", 
               requires = c("POUMM"), 
               sources = c(paste0(simName, ".R"), 
                           "../../defineSimReplication_ATSSeG0.R"),
               bsub.W = "24:00",
               bsub.mem = 6000,
               bsub.other = "-R beta ",
               sleep.every = 200, sleep.secs = 60)
}



execName <- "exec_20170428"

if(FALSE) {
  # transfer files to euler
  system('rsync -cv -r --progress ~/Documents/Bio/Projects/POUMM_1.3.0.tar.gz vmitov@euler:~/')
  #system('rsync -cv -r --progress ~/Documents/Bio/Projects/patherit_0.7.tar.gz vmitov@euler:~/')
  #system('rsync -cv -r --progress ~/Documents/Bio/Projects/benchtable_1.0.tar.gz vmitov@euler:~/')
  
  system(paste0('rsync -cv -r --progress ', '../defineSimReplication_ATSSeG0.R',  
                ' vmitov@euler:~/poummBayesianValidation/DATA2/'))
  
  #system('rsync -cv -r --progress vmitov@euler:~/POUMM_1.2.0_R_x86_64-slackware-linux-gnu.tar.gz .')
  
  
  system(paste0('rsync -cv -r --progress ', simName, '.R',  ' vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  system(paste0('rsync -cv -r --progress ', simName,  '.RData vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  system(paste0('rsync -cv -r --progress ', treeFileName,  ' vmitov@euler:~/poummBayesianValidation/DATA2/', simName, '/', execName, '/'))
  system(paste0('rsync -cv -r --progress j_* vmitov@euler:~/poummBayesianValidation/DATA2/', 
                simName, '/', execName, '/'))
  
} 

if(FALSE) {
  # bring only job files from euler 
  system(paste0('rsync -cv -r --remove-source-files --progress vmitov@euler:~/poummBayesianValidation/DATA/', simName, '/', execName, '/job* ', './', execName))
  
  # bring all files in simulation execution
  system(paste0('rsync -cv -r --remove-source-files --progress vmitov@euler:~/poummBayesianValidation/DATA/', simName, '/', execName, '/* ', './', execName))
}
