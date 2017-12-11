/*
 *  ThreePointMultivariate.h
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
#ifndef ParallelPruning_ThreePointMultivariate_H_
#define ParallelPruning_ThreePointMultivariate_H_

#include "splittree.h"
#include <armadillo>

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
namespace splittree {
template<class PruningTree>
class ThreePointMultivariate: public TraversalSpecification<PruningTree> {
public:
  typedef TraversalSpecification<PruningTree> BaseType;
  typedef PruningTree TreeType;
  typedef vec StateType;
  typedef vec ParameterType;

  // define fields as public in order to access them easily from R.
  arma::mat X, Y;
  arma::vec tTransf;
  arma::mat hat_mu_Y, tilde_mu_X_prime;
  arma::vec lnDetV, p, Q;

  ThreePointMultivariate(Rcpp::List const& tree): BaseType(tree) {
    this->tTransf = arma::vec(this->ref_tree_.num_nodes() - 1);
    this->lnDetV = arma::vec(this->ref_tree_.num_nodes(), arma::fill::zeros);
    this->p = arma::vec(this->ref_tree_.num_nodes(), arma::fill::zeros);
    this->Q = arma::vec(this->ref_tree_.num_nodes(), arma::fill::zeros);
  };

  void set_X_and_Y(arma::mat const& X, arma::mat const& Y) {
    if(X.n_rows != this->ref_tree_.num_tips() || Y.n_rows != this->ref_tree_.num_tips()) {
      Rcpp::stop("The matrices X and Y must have the same number of rows as V.");
    } else {
      this->X = X; this->Y = Y;

      this->hat_mu_Y = arma::mat(this->ref_tree_.num_nodes(), Y.n_cols, arma::fill::zeros);
      this->tilde_mu_X_prime = arma::mat(this->ref_tree_.num_nodes(), X.n_cols, arma::fill::zeros);
    }
  }

  StateType StateAtRoot() const {
    vec res(2);
    res[0] = this->lnDetV(this->ref_tree_.num_nodes()-1);
    res[1] = this->Q(this->ref_tree_.num_nodes()-1);
    return res;
  }

  inline void initSpecialData() {
    hat_mu_Y.fill(0);
    tilde_mu_X_prime.fill(0);
    lnDetV.fill(0);
    p.fill(0);
    Q.fill(0);
  }

  void InitNode(uint i) {
    tTransf[i] = this->ref_tree_.LengthOfBranch(i);
  }

  void VisitNode(uint i) {
    if(i < this->ref_tree_.num_nodes()) {
      // branch leading to a tip
      lnDetV[i] = log(tTransf[i]);
      p[i] = 1 / tTransf[i];
      hat_mu_Y.row(i) = Y.row(i);
      tilde_mu_X_prime.row(i) = X.row(i);
      Q[i] = arma::dot(X.row(i), Y.row(i)) / tTransf[i];
    } else {
      hat_mu_Y.row(i) /= p[i];
      tilde_mu_X_prime(i) /= p[i];
      Q[i] -= tTransf[i]*p[i]*p[i] / (1 + tTransf[i]*p[i]) *
        arma::dot(tilde_mu_X_prime.row(i), hat_mu_Y.row(i));
      lnDetV[i] += log(1 + tTransf[i]*p[i]);
      p[i] /= (1 + tTransf[i]*p[i]);
    }
  }

  void PruneNode(uint i, uint iParent) {
    hat_mu_Y.row(iParent) += p[i]*hat_mu_Y.row(i);
    tilde_mu_X_prime.row(iParent) += p[i]*tilde_mu_X_prime.row(i);
    lnDetV[iParent] += lnDetV[i];
    p[iParent] += p[i];
    Q[iParent] += Q[i];
  }
};
}
#endif // ParallelPruning_ThreePointMultivariate_H_
