/**
 *  AbcPMM.h
 *  SPLITT
 *
 * Copyright 2018 Venelin Mitov
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
#ifndef AbcPMM_H_
#define AbcPMM_H_

#include "./SPLITT.h"
#include "./NumericTraitData.h"
#include <iostream>
#include <cmath>

namespace PMMUsingSPLITT {

using namespace SPLITT;

template<class Tree>
class AbcPMM: public TraversalSpecification<Tree> {

public:
  typedef AbcPMM<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  double sigmae2, sigma2;
  vec x;
  vec a, b, c;

  AbcPMM(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.x_.size() != this->ref_tree_.num_tips()) {
      std::ostringstream oss;
      oss<<"The vector x must be the same length as the number of tips ("<<
        this->ref_tree_.num_tips()<<"), but were"<<input_data.x_.size()<<".";
      throw std::invalid_argument(oss.str());
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->x = At(input_data.x_, ordNodes);
      this->a = vec(this->ref_tree_.num_nodes());
      this->b = vec(this->ref_tree_.num_nodes());
      this->c = vec(this->ref_tree_.num_nodes());
    }
  };

  StateType StateAtRoot() const {
    vec res(3);
    res[0] = a[this->ref_tree_.num_nodes() - 1];
    res[1] = b[this->ref_tree_.num_nodes() - 1];
    res[2] = c[this->ref_tree_.num_nodes() - 1];
    return res;
  };

  void SetParameter(ParameterType const& par) {
    if(par.size() != 2) {
      throw std::invalid_argument(
      "The par vector should be of length 2 with \
      elements corresponding to sigma2 and sigmae2.");
    }
    if(par[0] <= 0 || par[1] <= 0) {
      throw std::logic_error("The parameters sigma2 and sigmae2 should be positive.");
    }
    this->sigma2 = par[0];
    this->sigmae2 = par[1];
  }

  inline void InitNode(uint i) {
    
    if(i < this->ref_tree_.num_tips()) {
      a[i] = -0.5 / sigmae2;  
      b[i] = x[i] / sigmae2;
      c[i] = -0.5 * (x[i]*b[i] + log(2*G_PI*sigmae2));
    } else {
      a[i] = b[i] = c[i] = 0;
    }
  }

  inline void VisitNode(uint i) {
    double t = this->ref_tree_.LengthOfBranch(i);
    
    double d = 1 - 2*a[i]*sigma2*t;
    // the order is important here because for c[i] we use the previous values 
    // of a[i] and b[i].
    c[i] = c[i] - 0.5*log(d) + 0.5*b[i]*b[i]*sigma2*t/d;
    a[i] /= d;
    b[i] /= d;
  }

  inline void PruneNode(uint i, uint j) {
    a[j] = a[j] + a[i];
    b[j] = b[j] + b[i];
    c[j] = c[j] + c[i];
  }

};
}

#endif //AbcPMM_H_
