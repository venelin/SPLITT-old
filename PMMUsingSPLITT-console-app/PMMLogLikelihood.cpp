/*
 *  PMMLogLikelihood.cpp
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

#include "SPLITT.h"
#include "./NumericTraitData.h"
#include "./AbcPMM.h"

#include <iostream>

using namespace SPLITT;
using namespace PMMUsingSPLITT;

typedef TraversalTask< AbcPMM<OrderedTree<std::string, double>> > ParallelPruningAbcPMM;

int main(int argc, char* argv[]) {
  
  // Will be false if the program is run without any arguments. Pass any argument
  // after the program name and verbose will be set to true.
  bool verbose = (argc > 1); 
  
  // will be using std::string, std::vector, std::cin and std::cout quite a lot.
  using namespace std;
  
  cout<<"Hello from the SPLITT PMM example!"<<endl;
 
  cout<<"Reading the input tree..."<<endl;
  
  // Read the number of nodes and tips
  uint M, N;
  cin>>M>>N;
  
  // read the tree
  vector<string> daughters(M-1);
  vector<string> parents(M-1);
  vec t(M-1);
  
  
  for(int i=0; i < M-1; i++) {
    cin>>daughters[i]>>parents[i]>>t[i];
    if(verbose) {
      cout<<daughters[i]<<" "<<parents[i]<<" "<<t[i]<<endl;
    }
  }
  
  cout<<"Reading the trait data..."<<endl;
  // read the trait data
  vector<string> tip_names(N);
  vec x(N);
  
  
  for(int i = 0; i < N; i++) {
    cin>>tip_names[i]>>x[i];
    if(verbose) {
      cout<<tip_names[i]<<" "<<x[i]<<endl;
    }
  }
  
  // create the data-object
  typename ParallelPruningAbcPMM::DataType data(tip_names, x);
  
  // Create the TraversalTask object: this will create the OrderedTree and the 
  // AbcPMM objects in the memory.
  ParallelPruningAbcPMM pruningTask(parents, daughters, t, data);
  
  // The tree objects reorders the nodes, so that the tips are indexed from
  // 0 to N-1, the internal nodes are from N to M-2, and the root is M-1.
  if(verbose) {
    cout<<"Node-name order in the OrderedTree:"<<endl;
    for(int i = 0; i < pruningTask.tree().num_nodes(); i++) {
      cout<<i<<". "<<pruningTask.tree().FindNodeWithId(i)<<endl;
    } 
  }
  
  // Check if OPENMP is enabled:
  std::cout<<"OpenMP-version: "<<pruningTask.algorithm().VersionOPENMP()<<std::endl;
  
  // model parameters
  double gM, sigma, sigmae;
  vec param(2); 
  
  
  cout<<"Main loop..."<<endl;
  
  // do this loop as long as parameters can be read from the standard input
  while( cin>>gM>>sigma>>sigmae ) {
    param[0] = sigma*sigma;
    param[1] = sigmae*sigmae;
  
    cout<<"  Calculating a, b, c at the root for sigma="<<sigma<<" and sigmae="<<sigmae<<std::endl;
    vec abc = pruningTask.TraverseTree(param, 0);
    cout<<"  a="<<abc[0]<<", b="<<abc[1]<<", c="<<abc[2]<<std::endl;
    
    cout<<"  LL(gM="<<gM<<", sigma="<<sigma<<", sigmae="<<sigmae<<"): "<<abc[0]*gM*gM + abc[1]*gM + abc[2]<<endl;
  } 
  
  cout<<"Main loop done."<<endl;
  // Exit politely
  std::cout<<"Good bye!"<<std::endl;
  return 0;
}


