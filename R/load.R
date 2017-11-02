#' @useDynLib ParallelPruning
NULL
# loading the ParallelPruning C++ modules
loadModule( "ParallelPruningTree", TRUE)
loadModule( "ParallelPruningTreeFindChildren", TRUE)
loadModule( "ParallelPruningTreeStringNodes", TRUE)
loadModule( "POUMM_abc", TRUE )
loadModule( "POUMM_lnDetV_Q_1d", TRUE )

