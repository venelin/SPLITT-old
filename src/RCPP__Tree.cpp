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

#include <R_ext/Rdynload.h>
#include <Rcpp.h>

#include "./SPLITT.h"
    
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

using namespace SPLITT;

SPLITT::Tree<uint, double>* CreateTree(Rcpp::List const& tree) {
  Rcpp::IntegerMatrix branches = tree["edge"];
  SPLITT::uvec br_0(branches.column(0).begin(), branches.column(0).end());
  SPLITT::uvec br_1(branches.column(1).begin(), branches.column(1).end());
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::Tree<uint, double>(br_0, br_1, t);
}

RCPP_MODULE(SPLITT__Tree) {
  Rcpp::class_< SPLITT::Tree<uint, double> > ( "SPLITT__Tree" )
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
  Rcpp::IntegerMatrix branches = tree["edge"];
  SPLITT::uvec br_0(branches.column(0).begin(), branches.column(0).end());
  SPLITT::uvec br_1(branches.column(1).begin(), branches.column(1).end());
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::OrderedTree<uint, double>(br_0, br_1, t);
}


RCPP_MODULE(SPLITT__OrderedTree) {
  Rcpp::class_< SPLITT::Tree<uint, double> > ( "SPLITT__Tree" )
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
  Rcpp::class_< SPLITT::OrderedTree<uint, double> >( "SPLITT__OrderedTree" )
    .derives< SPLITT::Tree<uint, double> > ( "SPLITT__Tree" )
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

RCPP_MODULE(SPLITT__OrderedTreeStringNodes) {
  Rcpp::class_< SPLITT::Tree<std::string, double> > ( "SPLITT__TreeStringNodes" )
  .property("num_nodes", &SPLITT::Tree<std::string, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<std::string, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<std::string, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<std::string, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<std::string, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<std::string, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<std::string, double>::FindChildren )
  .method("OrderNodes", &SPLITT::OrderedTree<std::string, double>::OrderNodes )
  ;
  Rcpp::class_< SPLITT::OrderedTree<std::string, double> >( "SPLITT__OrderedTreeStringNodes" )
    .derives< SPLITT::Tree<std::string, double> > ( "SPLITT__TreeStringNodes" )
    .factory<std::vector<std::string> const&, std::vector<std::string> const&,std::vector<double> const&>( &CreateOrderedTreeStringNodes )
    .property("num_levels", &SPLITT::OrderedTree<std::string, double>::num_levels )
    .property("num_parallel_ranges_prune", &SPLITT::OrderedTree<std::string, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &SPLITT::OrderedTree<std::string, double>::ranges_id_visit )
    .property("ranges_id_prune", &SPLITT::OrderedTree<std::string, double>::ranges_id_prune )
  ;
}
