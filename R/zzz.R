# zzz.R
# SPLITT
# 
# Copyright 2017 Venelin Mitov
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

## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}

# loading the RCPP C++ module
loadModule( "SPLITT__RCPP__BPM", TRUE )

loadModule( "SPLITT__RCPP__AbcPMM", TRUE )

loadModule( "SPLITT__RCPP__AbcPOUMM", TRUE )

loadModule( "SPLITT__RCPP__OrderedTreeStringNodes", TRUE )


# 
# 
# 
# library(testthat)
# library(SPLITT)
# library(microbenchmark)
# library(ape)
# library(POUMM)
# 
# lik_POUMM <- function(
#   poummObj,
#   g0, alpha, theta, sigma, sigmae, mode) {
#   
#   abc = poummObj$TraverseTree(c(alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae), mode)
#   
#   POUMM:::loglik_abc_g0_g0Prior(
#     abc = abc, g0Prior = NA,
#     g0 = g0, alpha = alpha, theta = theta, sigma = sigma)$loglik
# }
# 
# lik_POUMM_lnDetV_Q <- function(
#   poummObj,
#   g0, alpha, theta, sigma, sigmae, mode) {
#   
#   res <- poummObj$TraverseTree(c(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),mode)
#   
#   -1/2*(poummObj$tree$num_tips*log(2*pi) + 2*res[3] + res[1]+res[2])
# }
# 
# lik_POUMM_old <- function(
#   pruneInfo,
#   g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae) {
#   
#   POUMM::likPOUMMGivenTreeVTipsC(
#     integrator = pruneInfo$integrator,
#     g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae)
# }
# 
# context("PruneTree")
# set.seed(1)
# 
# EPS <- 1e-8
# N <- 100
# tree <- ape::rtree(N)
# 
# g0 <- 16
# alpha <- 2
# theta <- 4
# sigma <- .2
# sigmae <- .7
# se = rep(0, N); #rexp(N, 1/.01)
# z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)
# 
# # load("../../../poummBayesianValidation/DATA/Microbenchmark/Trees.RData")
# # tree<-trees$tree[[24]]
# # z<-trees$z[[24]]
# 
# 
# #ppt<-SPLITT:::ParallelPruningTree$new(tree)
# context("POUMM::pruneTree")
# pruneInfo <- POUMM::pruneTree(tree, z, se)
# 
# 
# context("SPLITT:::ParallelPruningThreePointPOUMM")
# poummlnDetVQ <- SPLITT:::ParallelPruningThreePointPOUMM$new(tree, z = z[1:N], se = se)
# 
# # test correct value
# test_that("POUMM lnDetV_Q", expect_lt(abs(
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21) -
#     lik_POUMM_old(pruneInfo,
#                   g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae
#     )), EPS))
# 
# 
# test_that("Equal parallel vs serial pruning", expect_equal(
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21),
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 1)))
# 
# 
# library(ape)
# library(SPLITT)
# library(microbenchmark)
# library(testthat)
# N <- 10000
# tree<-rtree(N)
# idx <- sample(which(tree$edge[, 2] > N), size=N/2)
# tree$edge.length[idx] <- 0
# #tree2 <- di2multi(tree)
# z <- rnorm(N)
# se <- rep(0, N)
# context("PruneTree")
# set.seed(1)
# 
# g0 <- 16
# alpha <- 2
# theta <- 4
# sigma <- .2
# sigmae <- .7
# 
# lik_POUMM_lnDetV_Q <- function(
#   poummObj,
#   g0, alpha, theta, sigma, sigmae, mode) {
#   
#   res <- poummObj$TraverseTree(c(g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae),mode)
#   
#   -1/2*(N*log(2*pi) + 2*res[3] + res[1]+res[2])
# }
# 
# poummlnDetVQ <- SPLITT:::ParallelPruningThreePointPOUMM$new(tree, z = z[1:N], se = se)
# 
# 
# 
# microbenchmark(
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 11),
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 12),
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 21),
#   lik_POUMM_lnDetV_Q(poummlnDetVQ,
#                      g0 = g0, alpha = alpha, theta = theta, sigma = sigma, sigmae = sigmae, mode = 22))
