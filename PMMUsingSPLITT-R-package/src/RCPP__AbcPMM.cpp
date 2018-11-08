/**
  *  RCPP__AbcPMM.cpp
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

#include <Rcpp.h>
#include "./AbcPMM.h"
    
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

using namespace SPLITT;
using namespace PMMUsingSPLITT;

typedef TraversalTask< AbcPMM<OrderedTree<uint, double>> > TraversalTaskAbcPMM;


TraversalTaskAbcPMM* CreateTraversalTaskAbcPMM(
    Rcpp::List const& tree, vec const& values) {
  
  Rcpp::IntegerMatrix branches = tree["edge"];
  uvec parents(branches.column(0).begin(), branches.column(0).end());
  uvec daughters(branches.column(1).begin(), branches.column(1).end());
  vec t = Rcpp::as<vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  uvec tip_names = Seq(uint(1), num_tips);
  
  typename TraversalTaskAbcPMM::DataType data(tip_names, values);
  
  return new TraversalTaskAbcPMM(parents, daughters, t, data);
}


// This will enable returning a copy of the `TraversalAlgorithm`-object stored in
// a `TraversalTaskAbcPMM` object to a R. This will be used in the MiniBenchmark
// R-function to check things like the OpenMP version used during compilation and
// the number of OpenMP threads at runtime. 
RCPP_EXPOSED_CLASS_NODECL(TraversalTaskAbcPMM::AlgorithmType)
  
RCPP_MODULE(PMMUsingSPLITT__TraversalTaskAbcPMM) {
  
  // Expose the properties VersionOPENMP and NumOmpThreads from the base 
  // TraversalAlgorithm class
  Rcpp::class_<TraversalTaskAbcPMM::AlgorithmType::ParentType> (
      "PMMUsingSPLITT__AbcPMM__TraversalAlgorithm"
    )
  .property( "VersionOPENMP",
             &TraversalTaskAbcPMM::AlgorithmType::ParentType::VersionOPENMP )
  .property( "NumOmpThreads",
             &TraversalTaskAbcPMM::AlgorithmType::ParentType::NumOmpThreads )
  ;

  // Expose the TraversalTaskAbcPMM::AlgorithmType specifying that it derives 
  // from the base TraversalAlgorithm class
  Rcpp::class_<TraversalTaskAbcPMM::AlgorithmType> (
      "PMMUsingSPLITT__AbcPMM__AlgorithmType"
    )
  .derives<TraversalTaskAbcPMM::AlgorithmType::ParentType>(
      "PMMUsingSPLITT__AbcPMM__TraversalAlgorithm"
    )
  ;
  
  // Finally, expose the TraversalTaskAbcPMM class - this is the main class in 
  // the module, which will be instantiated from R using the factory function
  // we've just written.
  Rcpp::class_<TraversalTaskAbcPMM>( "PMMUsingSPLITT__TraversalTaskAbcPMM" )
  // The <argument-type-list> MUST MATCH the arguments of the factory function 
  // defined above.
  .factory<Rcpp::List const&, vec const&>( &CreateTraversalTaskAbcPMM )
  // Expose the method that we will use to execute the TraversalTask
  .method( "TraverseTree", &TraversalTaskAbcPMM::TraverseTree )
  // Expose the algorithm property
  .property( "algorithm", &TraversalTaskAbcPMM::algorithm )
  ;
}

