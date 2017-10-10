#' @useDynLib ParallelPruning
NULL

#' Extract information for fast likelihood calculation using the
#'  breadth-first pruning algorithm.
#'
#' @param tree a phylo object
#'
#' @import Rcpp
#'
#' @return a list of objects used for likelihood evaluation
#'
#' @export
prepare <- function(tree) {

  edge <- tree$edge
  mode(edge) <- "integer"

  N <- length(tree$tip.label)                # number of tips
  M <- length(unique(as.vector(edge)))  # number of all nodes
  endingAt <- order(rbind(edge, c(0, N+1))[, 2])


  nonPrunedChildren <- rep(0, M)
  ee1 <- edge[, 1]
  while(length(ee1)) {
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    nonPrunedChildren[ee1[matchp]] <-
      nonPrunedChildren[ee1[matchp]] + 1
    ee1 <- ee1[-matchp]
  }

  tipsVector <- as.integer(c())
  tipsVectorIndex <- as.integer(c(0))

  edgeVector <- as.integer(c())
  edgeVectorIndex <- as.integer(c(0))

  # tips to be pruned; start by the tips in the initial tree
  tips <- as.integer(1:N)

  while(tips[1] != N+1) { # while the root has not become a tip itself
    tipsVectorIndex <-
      c(tipsVectorIndex,
        tipsVectorIndex[length(tipsVectorIndex)] + length(tips))

    # add the tips to be pruned to the tipsVector
    tipsVector <- c(tipsVector, tips)

    # indices in the edge-matrix of the to be pruned edges (pointing to tips)
    edgesToTips <- endingAt[tips]

    # empty the tips vector so it can be filled in with new tips to be pruned
    tips <- as.integer(c())

    edgesDone <- c()

    # vectorized update of the parent nodes.
    # resolving order of sibling edges being pruned at the same time
    while(length(edgesDone) != length(edgesToTips)) {
      if(length(edgesDone) == 0) {
        remainingParents <- unique(edge[edgesToTips, 1])
      } else {
        remainingParents <- unique(edge[edgesToTips[-edgesDone], 1])
      }

      # sib- edges that are first in edge[edgesToTips, 1] get served
      # first
      edgesNext <- match(remainingParents, edge[edgesToTips, 1])

      edgeVector <- c(edgeVector, edgesNext)

      # attach the index of the current last element of edgeVector
      # to edgeVectorIndex
      edgeVectorIndex <-
        c(edgeVectorIndex,
          edgeVectorIndex[length(edgeVectorIndex)] + length(edgesNext))

      # For the parent nodes, decrement the amount of non-visited
      # children
      nonPrunedChildren[edge[edgesToTips[edgesNext], 1]] <-
        nonPrunedChildren[edge[edgesToTips[edgesNext], 1]] - 1

      # those parent nodes with no more non-visited children will be
      # pruned next
      tips <-
        c(tips,
          edge[edgesToTips[edgesNext][nonPrunedChildren[edge[edgesToTips[edgesNext], 1]] == 0], 1])

      edgesToTips[edgesNext] <- NA
      edgesDone <- c(edgesDone, edgesNext)
    }
  }

  reordered <- reorderEdges(
    edge = edge, t = tree$edge.length,
    endingAt = endingAt,
    tipsVector = tipsVector, tipsVectorIndex = tipsVectorIndex,
    edgeVector = edgeVector, edgeVectorIndex = edgeVectorIndex,
    M = M, N = N, nLevels = length(tipsVectorIndex) - 1)

  c(list(M = M, N = N,
       tree = tree,
       endingAt = endingAt,
       tipsVector = tipsVector, tipsVectorIndex = tipsVectorIndex,
       nLevels = length(tipsVectorIndex) - 1,
       edgeVector = edgeVector, edgeVectorIndex = edgeVectorIndex),
    reordered)
}

reorderEdges <- function(edge, t,
                         endingAt,
                         tipsVector,
                         tipsVectorIndex,
                         edgeVector,
                         edgeVectorIndex,
                         M, N, nLevels) {

  parents <- edge[, 1]

  parents[parents == N+1] <- 2*M
  eReord <- parents
  tReord <- t
  reord <- ord <- 1:(length(t)+1)

  # edges pointing to tips
  edgesToTips <- endingAt[tipsVector[(tipsVectorIndex[1] + 1):tipsVectorIndex[2]]]

  # index in edgeVectorIndex
  jEVI <- 1

  lenEdgesDone <- 0
  while(lenEdgesDone != length(edgesToTips)) {
    edgesNext <- edgeVector[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]]

    edgeEnds <- edge[edgesToTips[edgesNext], 2]

    edgeEndsNew <- (edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]
    parents <- multiReplaceC(parents, edgeEnds, M+edgeEndsNew)

    eReord[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]] <-
      parents[edgesToTips[edgesNext]]
    eReord <- multiReplaceC(eReord, edgeEnds, M+edgeEndsNew)

    tReord[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]] <-
      t[edgesToTips[edgesNext]]

    reord[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]] <-
      ord[edgeEnds]

    jEVI <- jEVI + 1
    lenEdgesDone = lenEdgesDone + length(edgesNext)
  }

  # edges pointing to internal nodes which have become tips
  for(i in 2:nLevels) {
    edgesToTips <- endingAt[tipsVector[(tipsVectorIndex[i] + 1):tipsVectorIndex[i+1]]]

    lenEdgesDone <- 0
    while(lenEdgesDone != length(edgesToTips)) {
      edgesNext <- edgeVector[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]]
      #cat("edgesNext:"); print(edgesNext);

      edgeEnds <- edge[edgesToTips[edgesNext], 2]

      #cat("edgeEnds:"); print(edgeEnds)
      edgeEndsNew <- (edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]
      #cat("edgeEndsNew:"); print(edgeEndsNew)
      parents <- multiReplaceC(parents, edgeEnds, M+edgeEndsNew)

      eReord[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]] <- parents[edgesToTips[edgesNext]]
      eReord <- multiReplaceC(eReord, edgeEnds, M + edgeEndsNew)
      #cat("eReord:"); print(eReord)

      tReord[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]] <- t[edgesToTips[edgesNext]]
      #cat("tReord:"); print(tReord)

      reord[(edgeVectorIndex[jEVI] + 1):edgeVectorIndex[jEVI + 1]] <- ord[edgeEnds]
      #cat("reord:"); print(reord)

      jEVI <- jEVI + 1
      lenEdgesDone = lenEdgesDone + length(edgesNext)
    }
  }

  eReord <- eReord - M

  list(eReord = eReord, tReord = tReord, reord = reord)
}


# loading the ParallelPruning C++ modules
loadModule( "ParallelPruningAlgorithm", TRUE )
loadModule( "POUMM_Likelihood", TRUE )

