/*
 *  ThreePointUnivariate.h
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
#ifndef ParallelPruning_ThreePointUnivariate_H_
#define ParallelPruning_ThreePointUnivariate_H_

#include "./splittree.h"

namespace splittree {
// Calculate the |V| and quadratic quantities of the form Q=X'V^(-1)Y, for
// any covariance matrix V within the class of “3-point structured” matrices.
// X and Y must have the same number of rows as V but can have any number of
// columns. The algorithm also needs to calculate the following quantities
// p=1′V^(−1)1,  hat{mu}_Y=1′V^(−1)Y/p, and tilde{mu}_X'=X′V^(−1)1/p.
//
// The pruning procedure of the algorithm is described in the reference below.
// Here, we provide parallel pruning implementation of this algorithm.
//
// Reference: Lam Si Tung Ho and Cécile Ané. A Linear-Time Algorithm for
// Gaussian and Non-Gaussian Trait Evolution Models. SysBiol 2014.
template<class PruningTree>
class ThreePointUnivariate: public TraversalSpecification<PruningTree> {

public:
  typedef TraversalSpecification<PruningTree> BaseType;
  typedef PruningTree TreeType;
  typedef vec NodeStateType;
  typedef vec ParameterType;

  // public (unsafe) access to fields.
  vec X, Y;
  vec tTransf;
  vec hat_mu_Y, tilde_mu_X_prime;
  vec lnDetV, p, Q;

  ThreePointUnivariate(PruningTree const& tree): BaseType(tree) {
    this->tTransf = vec(this->ref_tree_.num_nodes() - 1);
    this->lnDetV = vec(this->ref_tree_.num_nodes(), 0);
    this->p = vec(this->ref_tree_.num_nodes(), 0);
    this->Q = vec(this->ref_tree_.num_nodes(), 0);
  };

  void set_X_and_Y(vec const& X, vec const& Y) {
    if(X.size() != this->ref_tree_.num_tips() || Y.size() != this->ref_tree_.num_tips()) {
      throw std::invalid_argument("The matrices X and Y must have the same number of rows as V.");
    } else {
      this->X = X; this->Y = Y;

      this->hat_mu_Y = vec(this->ref_tree_.num_nodes(), 0);
      this->tilde_mu_X_prime = vec(this->ref_tree_.num_nodes(), 0);
    }
  }

  NodeStateType StateAtRoot() const {
    vec res(2);
    res[0] = this->lnDetV[this->ref_tree_.num_nodes() - 1];
    res[1] = this->Q[this->ref_tree_.num_nodes() - 1];
    return res;
  }

  inline void InitNode(uint i) {
    hat_mu_Y[i] = tilde_mu_X_prime[i] = lnDetV[i] = p[i] = Q[i] = 0;
  }

  inline void VisitNode(uint i) {
    if(i < this->ref_tree_.num_tips()) {
      // branch leading to a tip
      lnDetV[i] = log(tTransf[i]);
      p[i] = 1 / tTransf[i];
      hat_mu_Y[i] = Y[i];
      tilde_mu_X_prime[i] = X[i];
      Q[i] = X[i] * Y[i] / tTransf[i];
    } else {
      hat_mu_Y[i] /= p[i];
      tilde_mu_X_prime[i] /= p[i];
      Q[i] -= tTransf[i]*p[i]*p[i] / (1 + tTransf[i]*p[i]) *
        tilde_mu_X_prime[i] * hat_mu_Y[i];
      lnDetV[i] += log(1 + tTransf[i]*p[i]);
      p[i] /= (1 + tTransf[i]*p[i]);
    }
  }

  inline void PruneNode(uint i, uint i_parent) {
    hat_mu_Y[i_parent] += p[i]*hat_mu_Y[i];
    tilde_mu_X_prime[i_parent] += p[i]*tilde_mu_X_prime[i];
    lnDetV[i_parent] += lnDetV[i];
    p[i_parent] += p[i];
    Q[i_parent] += Q[i];
  }
};
};

#endif // ParallelPruning_ThreePointUnivariate_H_
