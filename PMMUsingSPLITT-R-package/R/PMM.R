# PMM.R
# PMMUsingSPLITT
# 
# Copyright 2018 Venelin Mitov
# 
# This file is part of SPLITT: a generic C++ library for Serial and Parallel
# Lineage Traversal of Trees.
# 
# SPLITT is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# 
# SPLITT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with SPLITT.  If not, see
# <http://www.gnu.org/licenses/>.
# 
# @author Venelin Mitov


#' Calculate the PMM log-likelihood for a given tree, data and model parameters
#' @param x a numerical vector of size N, where N is the number of tips in tree
#' @param tree a phylo object
#' @param x0,sigma2,sigmae2 parameters of the PMM:
#' \describe{
#' \item{x0}{value at the root of tree excluding whte noise;}
#' \item{sigma2}{unit-time variance increment of the heritable component;}
#' \item{sigmae2}{variance of the non-heritable component.}
#' }
#' @param ord an integer vector denoting the pruning order. This should be the 
#' result from calling `reorder(tree, order = "postorder", index.only = TRUE)`, 
#' which is also set as a default value. Can be passed as argument to speed-up 
#' the calculation.
#' @return the log-likelihood value.
PMMLogLik <- function(
  x, tree, x0, sigma2, sigmae2, 
  ord = reorder(tree, order = "postorder", index.only = TRUE)) {
  
  # number of tips in the tree
  N <- length(tree$tip.label)
  # total number of nodes in the tree (tips, internal nodes and root node)
  M <- nrow(tree$edge) + 1L
  # state variables for each node
  a <- b <- c <- rep(0.0, M)

  # indices of the rows in tree$edge in pruning order
  if(is.null(ord)) {
    ord <- reorder(tree, order = "postorder", index.only = TRUE)
  }

  for(o in ord) {
    # daughter node
    i <- tree$edge[o, 2]
    # parent node
    j <- tree$edge[o, 1]
    # branch length
    t <- tree$edge.length[o]

    if(i <= N) {
      # initialize a tip node
      a[i] <- -0.5 / sigmae2
      b[i] <- x[i] / sigmae2
      c[i] <- -0.5 * (x[i]*b[i] + log(2*pi*sigmae2))
    }

    d = 1 - 2*a[i]*sigma2*t

    # the order is important
    c[i] <- c[i] - 0.5*log(d) + 0.5*b[i]*b[i]*sigma2*t/d
    a[i] <- a[i] / d
    b[i] <- b[i] / d
    
    
    a[j] <- a[j] + a[i]
    b[j] <- b[j] + b[i]
    c[j] <- c[j] + c[i]
  }

  # for phylo objects, N+1 denotes the root node
  a[N+1]*x0^2 + b[N+1]*x0 + c[N+1]
}

#' Calculate the PMM log-likelihood for a given tree, data and model parameters using the Rcpp module
#' @inheritParams PMMLogLik
#' @param cppObject a previously created object returned by \code{\link{NewPMMCppObject}}
#' @param mode an integer denoting the mode for traversing the tree, i.e. serial vs paralle.
#' 
#' @return the log-likelihood value.
PMMLogLikCpp <- function(x, tree, x0, sigma2, sigmae2, 
                         cppObject = NewPMMCppObject(x, tree),
                         mode = getOption("SPLITT.postorder.mode", 0)) {
  abc <- cppObject$TraverseTree(c(sigma2, sigmae2), mode)
  abc[1]*x0^2 + abc[2]*x0 + abc[3]
}


#' Create an instance of the Rcpp module for a given tree and trait data
#'
#' @inheritParams PMMLogLik
#' @return an object to be passed as argument of the \link{PMMLogLikCpp} function.
#' @seealso \link{PMMLogLikCpp}
NewPMMCppObject <- function(x, tree) {
  PMMUsingSPLITT__TraversalTaskAbcPMM$new(tree, x[1:length(tree$tip.label)])
}
