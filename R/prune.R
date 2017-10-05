#' @useDynLib ParallelPruning
NULL

#' Extract information for fast likelihood calculation using the
#'  breadth-first pruning algorithm.
#'
#' @param tree a phylo object
#' @param z Numeric vector with length(tree$tip.label) values
#' corresponding to tree$tip.label.
#' @param se Non-negative numerical or N-vector indicating known
#'  standard measurement error.
#'
#' @import Rcpp
#'
#' @return a list of objects used for likelihood evaluation
#'
#' @export
pruneTree <- function(tree, z, se = 0) {

  N <- length(tree$tip.label)                # number of tips
  M <- length(unique(as.vector(tree$edge)))  # number of all nodes
  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])

  edge <- tree$edge
  mode(tree$edge) <- "integer"

  if(length(se) != N) {
    se <- rep(se[1], N)
  }

  if(any(is.na(se) | se < 0)) {
    stop("All elements of se should be non-negative numbers.")
  }

  nonVisitedChildren <- rep(0, M)
  ee1 <- tree$edge[, 1]
  while(length(ee1)) {
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    nonVisitedChildren[ee1[matchp]] <-
      nonVisitedChildren[ee1[matchp]] + 1
    ee1 <- ee1[-matchp]
  }

  nodesVector <- as.integer(c())
  nodesIndex <- as.integer(c(0))

  unVector <- as.integer(c())
  unIndex <- as.integer(c(0))

  # start by pruning the tips
  nodes <- as.integer(1:N)

  while(nodes[1] != N+1) {
    nodesIndex <-
      c(nodesIndex, nodesIndex[length(nodesIndex)] + length(nodes))

    nodesVector <- c(nodesVector, nodes)

    # indices in the edge-matrix of the edges pointing to the to be
    # pruned nodes
    es <- endingAt[nodes]

    nodes <- as.integer(c())
    edgeEnds <- tree$edge[es, 2]


    # start with a un-vector that doesn't contain indices in es, so
    # es[-un] returns es intact.
    unAll <- length(es)+1
    lenUnAll <- 0L

    #update parent pifs
    while(lenUnAll != length(es)) {
      # indices in current es of the edges not yet added to
      # their parent node.
      un <- match(unique(tree$edge[es[-unAll], 1]), tree$edge[es, 1])

      unAll <- c(unAll, un)
      lenUnAll <- lenUnAll + length(un)

      # attach to the vector of such indices
      unVector <- c(unVector, un)

      # attach the index of the current last element of unVector
      # to unIndex
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))

      # For the parent nodes, decrement the amount of non-visited
      # children
      nonVisitedChildren[tree$edge[es[un], 1]] <-
        nonVisitedChildren[tree$edge[es[un], 1]] - 1

      # those parent nodes with no more non-visited children will be
      # pruned next
      nodes <-
        c(nodes,
          tree$edge[es[un][nonVisitedChildren[tree$edge[es[un], 1]] == 0], 1])

      es[un] <- NA
    }
  }

  reordered <- reorderEdges(
    edge = tree$edge, t = tree$edge.length,
    endingAt = endingAt, nodesVector = nodesVector,
    nodesIndex = nodesIndex, unVector = unVector, unIndex = unIndex,
    M = M, N = N, nLevels = length(nodesIndex) - 1)

  c(list(M = M, N = N,
       z = z, se = se, tree = tree,
       endingAt=endingAt,
       nodesVector=nodesVector, nodesIndex=nodesIndex,
       nLevels=length(nodesIndex)-1,
       unVector=unVector, unIndex=unIndex),
    reordered)
}

reorderEdges <- function(edge, t,
                         endingAt,
                         nodesVector,
                         nodesIndex,
                         unVector,
                         unIndex,
                         M, N, nLevels) {

  parents <- edge[, 1]

  parents[parents == N+1] <- 2*M
  eReord <- parents
  tReord <- t
  reord <- ord <- 1:N

  # edges pointing to tips
  es <- endingAt[nodesVector[(nodesIndex[1] + 1):nodesIndex[2]]]

  unJ <- 1

  lenUnAll <- 0
  while(lenUnAll != length(es)) {
    un <- unVector[(unIndex[unJ] + 1):unIndex[unJ + 1]]

    edgeEnds <- edge[es[un], 2]
    edgeEndsNew <- (unIndex[unJ] + 1):unIndex[unJ + 1]
    parents <- multiReplaceC(parents, edgeEnds, M+edgeEndsNew)

    eReord[(unIndex[unJ] + 1):unIndex[unJ + 1]] <- parents[es[un]]
    eReord <- multiReplaceC(eReord, edgeEnds, M+edgeEndsNew)
    tReord[(unIndex[unJ] + 1):unIndex[unJ + 1]] <- t[es[un]]

    reord[(unIndex[unJ] + 1):unIndex[unJ + 1]] <- ord[edgeEnds]

    unJ <- unJ + 1
    lenUnAll = lenUnAll + length(un)
  }

  # edges pointing to internal nodes
  for(i in 2:nLevels) {
    es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i+1]]]

    lenUnAll <- 0
    while(lenUnAll != length(es)) {
      un <- unVector[(unIndex[unJ] + 1):unIndex[unJ + 1]]

      edgeEnds <- edge[es[un], 2]

      edgeEndsNew <- (unIndex[unJ] + 1):unIndex[unJ + 1]
      parents <- multiReplaceC(parents, edgeEnds, M+edgeEndsNew)

      eReord[(unIndex[unJ] + 1):unIndex[unJ + 1]] <- parents[es[un]]
      eReord <- multiReplaceC(eReord, edgeEnds, M+edgeEndsNew)
      tReord[(unIndex[unJ] + 1):unIndex[unJ + 1]] <- t[es[un]]

      unJ <- unJ + 1
      lenUnAll = lenUnAll + length(un)
    }
  }

  eReord <- eReord - M

  list(eReord = eReord, tReord = tReord, reord = reord)
}


# loading the ImplementationCPP C++ module
loadModule( "ImplementationCPPParallelPruning", TRUE )

#' Extract information for fast likelihood calculation using the
#'  breadth-first pruning algorithm.
#'
#' @param tree a phylo object
#' @param z Numeric vector with length(tree$tip.label) values
#' corresponding to tree$tip.label.
#' @param se Non-negative numerical or N-vector indicating known
#'  standard measurement error.
#'
#' @import Rcpp
#'
#' @return a list of objects used for likelihood evaluation
#'
#' @export
pruneTree2 <- function(tree, z, se = 0) {

  N <- length(tree$tip.label)                # number of tips
  M <- length(unique(as.vector(tree$edge)))  # number of all nodes
  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])

  edge <- tree$edge
  mode(tree$edge) <- "integer"

  if(length(se) != N) {
    se <- rep(se[1], N)
  }

  if(any(is.na(se) | se < 0)) {
    stop("All elements of se should be non-negative numbers.")
  }

  nonVisitedChildren <- rep(0, M)
  ee1 <- tree$edge[, 1]
  while(length(ee1)) {
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    nonVisitedChildren[ee1[matchp]] <-
      nonVisitedChildren[ee1[matchp]] + 1
    ee1 <- ee1[-matchp]
  }

  nodesVector <- as.integer(c())
  nodesIndex <- as.integer(c(0))

  unVector <- as.integer(c())
  unIndex <- as.integer(c(0))

  # start by pruning the tips
  nodes <- as.integer(1:N)

  while(nodes[1] != N+1) {
    nodesIndex <-
      c(nodesIndex, nodesIndex[length(nodesIndex)] + length(nodes))

    nodesVector <- c(nodesVector, nodes)

    # indices in the edge-matrix of the edges pointing to the to be
    # pruned nodes
    es <- endingAt[nodes]

    nodes <- as.integer(c())
    edgeEnds <- tree$edge[es, 2]


    # start with a un-vector that doesn't contain indices in es, so
    # es[-un] returns es intact.
    unAll <- length(es)+1
    lenUnAll <- 0L

    #update parent pifs
    while(lenUnAll != length(es)) {
      # indices in current es of the edges not yet added to
      # their parent node.
      un <- match(unique(tree$edge[es[-unAll], 1]), tree$edge[es, 1])

      unAll <- c(unAll, un)
      lenUnAll <- lenUnAll + length(un)

      # attach to the vector of such indices
      unVector <- c(unVector, un)

      # attach the index of the current last element of unVector
      # to unIndex
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))

      # For the parent nodes, decrement the amount of non-visited
      # children
      nonVisitedChildren[tree$edge[es[un], 1]] <-
        nonVisitedChildren[tree$edge[es[un], 1]] - 1

      # those parent nodes with no more non-visited children will be
      # pruned next
      nodes <-
        c(nodes,
          tree$edge[es[un][nonVisitedChildren[tree$edge[es[un], 1]] == 0], 1])

      es[un] <- NA
    }
  }

  # create an ImplementationCPP object and initialize it with the pruning informaiton

  implementationCPP <- ImplementationCPP$new()
  implementationCPP$setPruningInfo(
    z, se, tree$edge, tree$edge.length,
    M, N, endingAt, nodesVector, nodesIndex, unVector, unIndex)

  list(M = M, N = N,
       z = z, se = se, tree = tree,
       endingAt=endingAt,
       nodesVector=nodesVector, nodesIndex=nodesIndex,
       nLevels=length(nodesIndex)-1,
       unVector=unVector, unIndex=unIndex,
       implementationCPP = implementationCPP)
}


#' For each i a[i] with b[i] in x
#' @note This R-implementation has a quadratic complexity
#' (length(x)*length(a)). This function has a faster implementation
#'  in Rcpp called multiReplaceC
#' @param x a numeric vector
#' @param a,b numeric vectors of the same length
#' @return The resulting vector x after the replacement
multiReplace <- function(x, a, b) {
  if(length(a) == 0) {
    stop("The vectors a and b should have at least 1 element")
  }
  if(length(a) != length(b)) {
    stop("The vectors a and b should have the same length")
  }
  for(i in 1:length(a)) {
    x[x==a[i]] <- b[i]
  }
  x
}

