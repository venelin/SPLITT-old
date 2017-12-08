#' @useDynLib SPLiTTree
NULL
# loading the SPLiTTree C++ modules
loadModule( "Tree1", TRUE)
loadModule( "OrderedTree1", TRUE)
loadModule( "OrderedTreeStringNodes", TRUE)
loadModule( "ParallelPruningAbcPOUMM", TRUE )
loadModule( "ParallelPruningThreePointPOUMM", TRUE )
