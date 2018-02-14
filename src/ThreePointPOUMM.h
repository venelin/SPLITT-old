/*
 *  ThreePointPOUMM.h
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

#ifndef ThreePointPOUMM_H_
#define ThreePointPOUMM_H_

#include "ThreePointUnivariate.h"
#include "NumericTraitData.h"

namespace splittree {

template<class Tree>
class ThreePointPOUMM: public ThreePointUnivariate<Tree> {

public:
  typedef ThreePointPOUMM<Tree> MyType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef ThreePointUnivariate<TreeType> BaseType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  // univariate trait vector
  splittree::vec z;
  splittree::vec h, u;
  double g0, alpha, alpha_x_2, theta, g0_theta, sigma, sigma2, sigma2_div_alpha_x_2, sigmae, sigmae2, e2alphaT, sum_u;
  double T; // tree height


  ThreePointPOUMM(
    TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.z_.size() != this->ref_tree_.num_tips() ||
       input_data.se_.size() != this->ref_tree_.num_tips()) {
      throw std::invalid_argument("ERR:01201:SPLiTTree:ThreePointPOUMM.h:ThreePointPOUMM:: The vectors z and se must be the same length as the number of tips.");
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->z = At(input_data.z_, ordNodes);
      vec X(this->ref_tree_.num_tips());
      this->set_X_and_Y(X, X);

      // A root-to-node distance vector in the order of pruning processing
      h.resize(this->ref_tree_.num_nodes() - 1);
      for(int i = this->ref_tree_.num_nodes() - 2; i >= 0; i--) {
        h[i] = h[this->ref_tree_.FindIdOfParent(i)] + this->ref_tree_.LengthOfBranch(i);
      }

      this->T = *std::max_element(h.begin(), h.begin()+this->ref_tree_.num_tips());
      this->u = splittree::vec(this->ref_tree_.num_tips());
      for(int i = 0; i < this->ref_tree_.num_tips(); i++) u[i] = T - h[i];
      sum_u = 0;
      for(auto uu : u) sum_u += uu;
    }
  }

  void SetParameter(ParameterType const& par) {
    if(par.size() != 5) {
      throw std::invalid_argument(
          "ERR:01211:SPLiTTree:ThreePointPOUMM.h:SetParameter:: The par vector should be of length 5 with \
      elements corresponding to g0, alpha, theta, sigma and sigmae.");
    }
    if(par[1] < 0 || par[3] < 0 || par[4] < 0) {
      throw std::logic_error("ERR:01212:SPLiTTree:ThreePointPOUMM.h:SetParameter:: The parameters alpha, sigma and sigmae should be non-negative.");
    }
    this->g0 = par[0];
    this->alpha = par[1];
    this->alpha_x_2 = 2*alpha;
    this->theta = par[2];
    this->sigma = par[3];
    this->sigma2 = sigma*sigma;
    this->sigmae = par[4];
    this->sigmae2 = sigmae*sigmae;
    this->sigma2_div_alpha_x_2 = sigma2/alpha_x_2;
    this->g0_theta = g0 - theta;
    this->e2alphaT = exp(-2*alpha*T);
  }

  StateType StateAtRoot() const {
    vec res = BaseType::StateAtRoot();
    res.push_back(alpha*sum_u);
    return res;
  }

  inline void InitNode(uint i) {
    ThreePointUnivariate<TreeType>::InitNode(i);
    if(i < this->ref_tree_.num_nodes() - 1) {
      // if an internal node or a tip, transform the branch length leading to this tip
      uint iParent = this->ref_tree_.FindIdOfParent(i);
      // tTransf[i] =
      //   sigma*sigma/(2*alpha) * ( (1-exp(-2*alpha*h[i]))*exp(-2*alpha*(T-h[i])) -
      //     (1-exp(-2*alpha*h[iParent]))*exp(-2*alpha*(T-h[iParent])) );

      double ealphahi = exp(alpha*h[i]);

      this->tTransf[i] = sigma2_div_alpha_x_2 *
        (e2alphaT*(ealphahi*ealphahi - exp(alpha_x_2*h[iParent])));

      if(i < this->ref_tree_.num_tips()) {
        //double mu = exp(-alpha*h[i])*g0 + (1-exp(-alpha*h[i]))*theta;
        // double mu = g0 / ealphahi + (1 - 1/ealphahi)*theta;
        double mu = theta + g0_theta/ealphahi;
        double ealphaui = exp(-alpha*u[i]);
        this->X[i] = this->Y[i] = (z[i] - mu)*ealphaui;
        this->tTransf[i] += sigmae2 * ealphaui*ealphaui;
      }
    }
  }
};

typedef TraversalTask<
  ThreePointPOUMM<OrderedTree<uint, double>> > ParallelPruningThreePointPOUMM;

}

#endif // ThreePointPOUMM_H_
