#' @useDynLib ParallelPruning
NULL
# loading the ParallelPruning C++ modules
loadModule( "Tree1", TRUE)
loadModule( "ParallelPruningTree1", TRUE)
loadModule( "ParallelPruningTreeStringNodes", TRUE)
loadModule( "ParallelPruningTreeFindChildren", TRUE)
loadModule( "ParallelPruningAbcPOUMM", TRUE )
loadModule( "ParallelPruningThreePointPOUMM", TRUE )

loadModule( "X", TRUE)
loadModule( "Y" , TRUE)
