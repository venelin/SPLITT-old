/*
 *  ppaRcpp.h
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
#include<iostream>
#include "ThreePointPOUMM.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// BEGIN: Needed for r-devel (R 3.4)
void R_init_SPLITT(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_SPLITT(DllInfo *info) {
  /* Release resources. */
}
// END Needed for r-devel (R 3.4)

SPLITT::Tree<uint, double>* CreateTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::Tree<uint, double>(br_0, br_1, t);
}

RCPP_MODULE(Tree1) {
  Rcpp::class_<SPLITT::Tree<uint, double>> ( "Tree1" )
  .factory<Rcpp::List const&>( &CreateTree )
  .property("num_nodes", &SPLITT::Tree<uint, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<uint, double>::FindChildren )
  ;
}

SPLITT::OrderedTree<uint, double>* CreateOrderedTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::OrderedTree<uint, double>(br_0, br_1, t);
}


RCPP_MODULE(OrderedTree1) {
  Rcpp::class_<SPLITT::Tree<uint, double>> ( "Tree1" )
  .factory<Rcpp::List const&>( &CreateTree )
  .property("num_nodes", &SPLITT::Tree<uint, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &SPLITT::Tree<uint, double>::OrderNodes )
  ;
  Rcpp::class_<SPLITT::OrderedTree<uint, double>>( "OrderedTree1" )
    .derives<SPLITT::Tree<uint, double>> ( "Tree1" )
    .factory<Rcpp::List const&>( &CreateOrderedTree )
    .property("num_levels", &SPLITT::OrderedTree<uint, double>::num_levels )
    .property("num_parallel_ranges_prune", &SPLITT::OrderedTree<uint, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &SPLITT::OrderedTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &SPLITT::OrderedTree<uint, double>::ranges_id_prune )
  ;
}

SPLITT::OrderedTree<std::string, double>* CreateOrderedTreeStringNodes(
    std::vector<std::string> const& br_0,
    std::vector<std::string> const& br_1,
    std::vector<double> const& t) {

  return new SPLITT::OrderedTree<std::string, double>(br_0, br_1, t);
}

RCPP_MODULE(OrderedTreeStringNodes) {
  Rcpp::class_<SPLITT::Tree<std::string, double>> ( "TreeStringNodes" )
  .property("num_nodes", &SPLITT::Tree<std::string, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<std::string, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<std::string, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<std::string, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<std::string, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<std::string, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<std::string, double>::FindChildren )
  .method("OrderNodes", &SPLITT::OrderedTree<std::string, double>::OrderNodes )
  ;
  Rcpp::class_<SPLITT::OrderedTree<std::string, double>>( "OrderedTreeStringNodes" )
    .derives<SPLITT::Tree<std::string, double>> ( "TreeStringNodes" )
    .factory<std::vector<std::string> const&, std::vector<std::string> const&,std::vector<double> const&>( &CreateOrderedTreeStringNodes )
    .property("num_levels", &SPLITT::OrderedTree<std::string, double>::num_levels )
    .property("num_parallel_ranges_prune", &SPLITT::OrderedTree<std::string, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &SPLITT::OrderedTree<std::string, double>::ranges_id_visit )
    .property("ranges_id_prune", &SPLITT::OrderedTree<std::string, double>::ranges_id_prune )
  ;
}

SPLITT::ParallelPruningThreePointPOUMM* CreateParallelPruningThreePointPOUMM(
    Rcpp::List const& tree, SPLITT::vec const& z, SPLITT::vec const& se) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec node_names = SPLITT::Seq(uint(1), num_tips);
  typename SPLITT::ParallelPruningThreePointPOUMM::DataType data(node_names, z, se);
  return new SPLITT::ParallelPruningThreePointPOUMM(br_0, br_1, t, data);
}

RCPP_EXPOSED_CLASS_NODECL(SPLITT::ParallelPruningThreePointPOUMM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType)


RCPP_MODULE(ParallelPruningThreePointPOUMM) {
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree> ( "Tree" )
  .property("num_nodes", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::num_nodes )
  .property("num_tips", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::num_tips )
  .method("LengthOfBranch", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM::TreeType>( "OrderedTree" )
    .derives<SPLITT::ParallelPruningThreePointPOUMM::TreeType::Tree> ( "Tree" )
    .method("RangeIdPruneNode", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::RangeIdVisitNode )
    .property("num_levels", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::num_levels )
    .property("ranges_id_visit", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &SPLITT::ParallelPruningThreePointPOUMM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType> ( "TraversalAlgorithm" )
    .property( "VersionOPENMP", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType::num_threads )
  ;    
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType> ( "PostOrderTraversal" )
    .derives<SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType> ( "TraversalAlgorithm" )
    .method( "ModeAutoStep", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &SPLITT::ParallelPruningThreePointPOUMM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType>( "ThreePointUnivariate" )
  .field( "Q", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::Q )
  .field( "lnDetV_all", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::lnDetV )
  .field( "p", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::p )
  .field( "X", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::X )
  .field( "Y", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::Y )
  .field( "tTransf", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::tTransf )
  .field( "hat_mu_Y", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::hat_mu_Y )
  .field( "tilde_mu_X_prime", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::tilde_mu_X_prime )
  ;
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType>( "ThreePointPOUMM" )
    .derives<SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType>("ThreePointUnivariate")
    .field( "h", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::h )
    .field( "u", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::u )
    .field( "alpha", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::alpha )
    .field( "sigma", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::sigma )
    .field( "sigmae", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::sigmae )
    .field( "g0", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::g0 )
    .field( "theta", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::theta )
    .field( "T", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::T )
    .field( "e2alphaT", &SPLITT::ParallelPruningThreePointPOUMM::TraversalSpecificationType::e2alphaT )
  ;
  Rcpp::class_<SPLITT::ParallelPruningThreePointPOUMM>( "ParallelPruningThreePointPOUMM" )
    .factory<Rcpp::List const&, SPLITT::vec const&, SPLITT::vec const&>(&CreateParallelPruningThreePointPOUMM)
    .method( "TraverseTree", &SPLITT::ParallelPruningThreePointPOUMM::TraverseTree )
    .property( "tree", &SPLITT::ParallelPruningThreePointPOUMM::tree )
    .property( "spec", &SPLITT::ParallelPruningThreePointPOUMM::spec )
    .property( "algorithm", &SPLITT::ParallelPruningThreePointPOUMM::algorithm )
  ;
}