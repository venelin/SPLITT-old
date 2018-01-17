/*
 *  ppaRcpp.h
 *  SPLiTTree
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of SPLiTTree: a generic C++ library for Serial and Parallel
 * Lineage Traversal of Trees.
 *
 * SPLiTTree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * SPLiTTree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SPLiTTree.  If not, see
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
void R_init_SPLiTTree(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_SPLiTTree(DllInfo *info) {
  /* Release resources. */
}
// END Needed for r-devel (R 3.4)

splittree::Tree<uint, double>* CreateTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  return new splittree::Tree<uint, double>(br_0, br_1, t);
}

RCPP_MODULE(Tree1) {
  Rcpp::class_<splittree::Tree<uint, double>> ( "Tree1" )
  .factory<Rcpp::List const&>( &CreateTree )
  .property("num_nodes", &splittree::Tree<uint, double>::num_nodes )
  .property("num_tips", &splittree::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &splittree::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &splittree::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &splittree::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &splittree::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &splittree::Tree<uint, double>::FindChildren )
  ;
}

splittree::OrderedTree<uint, double>* CreateOrderedTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  return new splittree::OrderedTree<uint, double>(br_0, br_1, t);
}


RCPP_MODULE(OrderedTree1) {
  Rcpp::class_<splittree::Tree<uint, double>> ( "Tree1" )
  .factory<Rcpp::List const&>( &CreateTree )
  .property("num_nodes", &splittree::Tree<uint, double>::num_nodes )
  .property("num_tips", &splittree::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &splittree::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &splittree::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &splittree::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &splittree::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &splittree::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &splittree::Tree<uint, double>::OrderNodes )
  ;
  Rcpp::class_<splittree::OrderedTree<uint, double>>( "OrderedTree1" )
    .derives<splittree::Tree<uint, double>> ( "Tree1" )
    .factory<Rcpp::List const&>( &CreateOrderedTree )
    .property("num_levels", &splittree::OrderedTree<uint, double>::num_levels )
    .property("num_parallel_ranges_prune", &splittree::OrderedTree<uint, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &splittree::OrderedTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &splittree::OrderedTree<uint, double>::ranges_id_prune )
  ;
}

splittree::OrderedTree<std::string, double>* CreateOrderedTreeStringNodes(
    std::vector<std::string> const& br_0,
    std::vector<std::string> const& br_1,
    std::vector<double> const& t) {

  return new splittree::OrderedTree<std::string, double>(br_0, br_1, t);
}

RCPP_MODULE(OrderedTreeStringNodes) {
  Rcpp::class_<splittree::Tree<std::string, double>> ( "TreeStringNodes" )
  .property("num_nodes", &splittree::Tree<std::string, double>::num_nodes )
  .property("num_tips", &splittree::Tree<std::string, double>::num_tips )
  .method("LengthOfBranch", &splittree::Tree<std::string, double>::LengthOfBranch )
  .method("FindNodeWithId", &splittree::Tree<std::string, double>::FindNodeWithId )
  .method("FindIdOfNode", &splittree::Tree<std::string, double>::FindIdOfNode )
  .method("FindIdOfParent", &splittree::Tree<std::string, double>::FindIdOfParent )
  .method( "FindChildren", &splittree::Tree<std::string, double>::FindChildren )
  .method("OrderNodes", &splittree::OrderedTree<std::string, double>::OrderNodes )
  ;
  Rcpp::class_<splittree::OrderedTree<std::string, double>>( "OrderedTreeStringNodes" )
    .derives<splittree::Tree<std::string, double>> ( "TreeStringNodes" )
    .factory<std::vector<std::string> const&, std::vector<std::string> const&,std::vector<double> const&>( &CreateOrderedTreeStringNodes )
    .property("num_levels", &splittree::OrderedTree<std::string, double>::num_levels )
    .property("num_parallel_ranges_prune", &splittree::OrderedTree<std::string, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &splittree::OrderedTree<std::string, double>::ranges_id_visit )
    .property("ranges_id_prune", &splittree::OrderedTree<std::string, double>::ranges_id_prune )
  ;
}

splittree::ParallelPruningThreePointPOUMM* CreateParallelPruningThreePointPOUMM(
    Rcpp::List const& tree, splittree::vec const& z, splittree::vec const& se) {
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  splittree::uvec node_names = splittree::Seq(uint(1), num_tips);
  typename splittree::ParallelPruningThreePointPOUMM::DataType data(node_names, z, se);
  return new splittree::ParallelPruningThreePointPOUMM(br_0, br_1, t, data);
}

RCPP_EXPOSED_CLASS_NODECL(splittree::ParallelPruningThreePointPOUMM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(splittree::ParallelPruningThreePointPOUMM::AlgorithmType)


RCPP_MODULE(ParallelPruningThreePointPOUMM) {
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM::TreeType::Tree> ( "Tree" )
  .property("num_nodes", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::num_nodes )
  .property("num_tips", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::num_tips )
  .method("LengthOfBranch", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::LengthOfBranch )
  .method("FindNodeWithId", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &splittree::ParallelPruningThreePointPOUMM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM::TreeType>( "OrderedTree" )
    .derives<splittree::ParallelPruningThreePointPOUMM::TreeType::Tree> ( "Tree" )
    .method("RangeIdPruneNode", &splittree::ParallelPruningThreePointPOUMM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &splittree::ParallelPruningThreePointPOUMM::TreeType::RangeIdVisitNode )
    .property("num_levels", &splittree::ParallelPruningThreePointPOUMM::TreeType::num_levels )
    .property("ranges_id_visit", &splittree::ParallelPruningThreePointPOUMM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &splittree::ParallelPruningThreePointPOUMM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType> ( "TraversalAlgorithm" )
    .property( "VersionOPENMP", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType::num_threads )
  ;    
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM::AlgorithmType> ( "PostOrderTraversal" )
    .derives<splittree::ParallelPruningThreePointPOUMM::AlgorithmType::ParentType> ( "TraversalAlgorithm" )
    .method( "ModeAutoStep", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &splittree::ParallelPruningThreePointPOUMM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType>( "ThreePointUnivariate" )
  .field( "Q", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::Q )
  .field( "lnDetV_all", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::lnDetV )
  .field( "p", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::p )
  .field( "X", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::X )
  .field( "Y", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::Y )
  .field( "tTransf", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::tTransf )
  .field( "hat_mu_Y", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::hat_mu_Y )
  .field( "tilde_mu_X_prime", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType::tilde_mu_X_prime )
  ;
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType>( "ThreePointPOUMM" )
    .derives<splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::BaseType>("ThreePointUnivariate")
    .field( "h", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::h )
    .field( "u", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::u )
    .field( "alpha", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::alpha )
    .field( "sigma", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::sigma )
    .field( "sigmae", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::sigmae )
    .field( "g0", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::g0 )
    .field( "theta", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::theta )
    .field( "T", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::T )
    .field( "e2alphaT", &splittree::ParallelPruningThreePointPOUMM::TraversalSpecificationType::e2alphaT )
  ;
  Rcpp::class_<splittree::ParallelPruningThreePointPOUMM>( "ParallelPruningThreePointPOUMM" )
    .factory<Rcpp::List const&, splittree::vec const&, splittree::vec const&>(&CreateParallelPruningThreePointPOUMM)
    .method( "TraverseTree", &splittree::ParallelPruningThreePointPOUMM::TraverseTree )
    .property( "tree", &splittree::ParallelPruningThreePointPOUMM::tree )
    .property( "spec", &splittree::ParallelPruningThreePointPOUMM::spec )
    .property( "algorithm", &splittree::ParallelPruningThreePointPOUMM::algorithm )
  ;
}