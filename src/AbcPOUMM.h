/*
 *  AbcPOUMM.h
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

#ifndef ABC_POUMM_H_
#define ABC_POUMM_H_

#include "./splittree.h"
#include "./NumericTraitData.h"

namespace splittree {

template<class Tree>
class AbcPOUMM: public TraversalSpecification<Tree> {

public:
  typedef AbcPOUMM<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> InputDataType;
  typedef vec NodeStateType;

  double alpha, theta, sigmae2, sigma2;
  vec z, se;
  vec a, b, c;

  AbcPOUMM(TreeType const& tree, InputDataType const& input_data):
    BaseType(tree) {

    if(input_data.z_.size() != this->ref_tree_.num_tips() ||
       input_data.se_.size() != this->ref_tree_.num_tips()) {
      throw std::invalid_argument("The vectors z and se must be the same length as the number of tips.");
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->z = At(input_data.z_, ordNodes);
      this->se = At(input_data.se_, ordNodes);
      this->a = vec(this->ref_tree_.num_nodes());
      this->b = vec(this->ref_tree_.num_nodes());
      this->c = vec(this->ref_tree_.num_nodes());
    }
  };

  NodeStateType StateAtRoot() const {
    vec res(3);
    res[0] = a[this->ref_tree_.num_nodes() - 1];
    res[1] = b[this->ref_tree_.num_nodes() - 1];
    res[2] = c[this->ref_tree_.num_nodes() - 1];
    return res;
  };

  void SetParameter(ParameterType const& par) {
    if(par.size() != 4) {
      throw std::invalid_argument(
      "The par vector should be of length 4 with \
      elements corresponding to alpha, theta, sigma and sigmae.");
    }
    if(par[0] < 0 || par[2] < 0 || par[3] < 0) {
      throw std::logic_error("The parameters alpha, sigma and sigmae should be non-negative.");
    }
    this->alpha = par[0];
    this->theta = par[1];
    this->sigma2 = par[2]*par[2];
    this->sigmae2 = par[3]*par[3];
  }

  inline void InitNode(uint i) {
    if(i < this->ref_tree_.num_tips()) {
      double sum_se2_sigmae2 = sigmae2 + se[i]*se[i];
      double z1 = z[i] - theta;
      a[i] = -0.5 / sum_se2_sigmae2;
      b[i] = z1 / sum_se2_sigmae2;
      c[i] = -0.5 * (M_LN_2PI  + z1 * b[i] + log(sum_se2_sigmae2));
    } else {
      a[i] = b[i] = c[i] = 0;
    }
  }

  inline void VisitNode(uint i) {

    double t = this->ref_tree_.LengthOfBranch(i);
    double talpha = t * alpha;
    double etalpha = exp(talpha);
    double e2talpha = etalpha * etalpha;
    double fe2talpha;
    if(alpha != 0) {
      fe2talpha = alpha / (1 - e2talpha);
    } else {
      fe2talpha = -0.5 / t;
    }
    double gutalphasigma2 = e2talpha + (a[i] * sigma2) / fe2talpha;

    c[i] = -0.5 * log(gutalphasigma2) - 0.25 * sigma2 * b[i] * b[i] /
      (fe2talpha - alpha + a[i] * sigma2) + talpha + c[i];
    b[i] = (etalpha * b[i]) / gutalphasigma2;
    a[i] /= gutalphasigma2;
  }

  inline void PruneNode(uint i, uint i_parent) {
    a[i_parent] += a[i];
    b[i_parent] += b[i];
    c[i_parent] += c[i];
  }

};

typedef TraversalTask<
  AbcPOUMM<OrderedTree<uint, double>> > ParallelPruningAbcPOUMM;
}
#endif //ABC_POUMM_H_
