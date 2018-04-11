#' @useDynLib SPLITT
NULL
# loading the SPLITT C++ modules
loadModule( "Tree1", TRUE)
loadModule( "OrderedTree1", TRUE)
loadModule( "OrderedTreeStringNodes", TRUE)
loadModule( "ParallelPruningThreePointPOUMM", TRUE )
