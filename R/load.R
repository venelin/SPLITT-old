#' @useDynLib ParallelPruning
NULL
# loading the ParallelPruning C++ modules
loadModule( "Tree1", TRUE)
loadModule( "PruningTree1", TRUE)
loadModule( "PruningTreeStringNodes", TRUE)
loadModule( "ParallelPruningAbcPOUMM", TRUE )
loadModule( "ParallelPruningThreePointPOUMM", TRUE )

loadModule( "X", TRUE)
loadModule( "Y" , TRUE)
