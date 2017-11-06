#' @useDynLib ParallelPruning
NULL
# loading the ParallelPruning C++ modules
loadModule( "Tree", TRUE)
loadModule( "ParallelPruningTree", TRUE)
loadModule( "ParallelPruningTreeStringNodes", TRUE)
loadModule( "ParallelPruningTreeFindChildren", TRUE)
loadModule( "ParallelPruningAbcPOUMM", TRUE )
loadModule( "ParallelPruningThreePointPOUMM", TRUE )

loadModule( "X", TRUE)
loadModule( "Y" , TRUE)
