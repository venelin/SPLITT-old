require(benchtable)
 require(POUMM)
simType <- "POUMM"
source('Ultrametric_N4000_ATSSeG0.R')
 source('../../defineSimReplication_ATSSeG0.R')

#user code 

# end of user code
args <- commandArgs(trailingOnly = TRUE)
f <- as.character(args[1])
table.file <- as.character(args[2])
table.name <- as.character(args[3])
ids <- as.integer(args[-(1:3)])
b <- doBenchJob(fname="replication", f=f, table.file=table.file, table.name=table.name, ids=ids)

