/**
  *  Rcpp.cpp
  *  SPLITT
  *
  * Copyright 2017 Venelin Mitov
  *
  * This file is part of SPLITT: a generic C++ library for Serial and Parallel
  * Lineage Traversal of Trees.
  *
  * SPLITT is free software: you can redistribute it and/or modify
  * it under the terms of the GNU Lesser General Public License as
  * published by the Free Software Foundation, either version 3 of
  * the License, or (at your option) any later version.
  *
  * SPLITT is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU Lesser General Public License for more details.
  *
  * You should have received a copy of the GNU Lesser General Public
  * License along with SPLITT.  If not, see
  * <http://www.gnu.org/licenses/>.
  *
  * @author Venelin Mitov
  */

#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
    
#include "./AbcPMM.h"
    
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace SPLITT;
using namespace PMMUsingSPLITT;

typedef TraversalTask< AbcPMM<OrderedTree<uint, double>> > ParallelPruningAbcPMM;


ParallelPruningAbcPMM* CreateParallelPruningAbcPMM(
    Rcpp::List const& tree, vec const& x) {
  arma::umat branches = tree["edge"];
  uvec br_0 = arma::conv_to<uvec>::from(branches.col(0));
  uvec br_1 = arma::conv_to<uvec>::from(branches.col(1));
  vec t = Rcpp::as<vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  uvec node_names = Seq(uint(1), num_tips);
  typename ParallelPruningAbcPMM::DataType data(node_names, x);
  return new ParallelPruningAbcPMM(br_0, br_1, t, data);
}

RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPMM::AlgorithmType)
  
RCPP_MODULE(PMMUsingSPLITT__AbcPMM) {
  Rcpp::class_<ParallelPruningAbcPMM::AlgorithmType::ParentType> ( "PMMUsingSPLITT__AbcPMM__TraversalAlgorithm" )
    .property( "VersionOPENMP", &ParallelPruningAbcPMM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &ParallelPruningAbcPMM::AlgorithmType::ParentType::NumOmpThreads )
  ;
  Rcpp::class_<ParallelPruningAbcPMM::AlgorithmType> ( "PMMUsingSPLITT__AbcPMM__AlgorithmType" )
    .derives<ParallelPruningAbcPMM::AlgorithmType::ParentType>( "PMMUsingSPLITT__AbcPMM__TraversalAlgorithm" )
  ;
  Rcpp::class_<ParallelPruningAbcPMM>( "PMMUsingSPLITT__AbcPMM" )
  .factory<Rcpp::List const&, vec const&>(&CreateParallelPruningAbcPMM)
  .method( "TraverseTree", &ParallelPruningAbcPMM::TraverseTree )
  .property( "algorithm", &ParallelPruningAbcPMM::algorithm )
  ;
}

