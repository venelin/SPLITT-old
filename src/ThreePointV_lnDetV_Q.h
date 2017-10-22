//   Copyright 2017 Venelin Mitov
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
//   limitations under the License.

#ifndef ParallelPruning_ThreePointV_lnDetV_Q_H_
#define ParallelPruning_ThreePointV_lnDetV_Q_H_

#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include "./ParallelPruningAlgorithm.h"
#include "./ppaRcpp.h"
#include <iostream>


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
namespace ppa {
class ThreePointV_lnDetV_Q {
public:
  typedef ParallelPruningMeta<ThreePointV_lnDetV_Q> MetaType;
  const MetaType meta;
  // define fields as public in order to access them easily from R.
  arma::mat X, Y;
  arma::vec tTransf;
  arma::mat hat_mu_Y, tilde_mu_X_prime;
  arma::vec lnDetV, p, Q;

  ThreePointV_lnDetV_Q(Rcpp::List const& tree):
    meta(ppaTreeFromRList(tree)) {

    this->tTransf = arma::vec(meta.M - 1);
    this->lnDetV = arma::vec(meta.M, arma::fill::zeros);
    this->p = arma::vec(meta.M, arma::fill::zeros);
    this->Q = arma::vec(meta.M, arma::fill::zeros);
  };

  void set_X_and_Y(arma::mat const& X, arma::mat const& Y) {
    if(X.n_rows != meta.N || Y.n_rows != meta.N) {
      Rcpp::stop("The matrices X and Y must have the same number of rows as V.");
    } else {
      this->X = X; this->Y = Y;

      this->hat_mu_Y = arma::mat(meta.M, Y.n_cols, arma::fill::zeros);
      this->tilde_mu_X_prime = arma::mat(meta.M, X.n_cols, arma::fill::zeros);
    }
  }

  double get_Q() const {
    return this->Q(meta.M-1);
  }

  double get_lnDetV() const {
    return this->lnDetV(meta.M-1);
  }

  inline void initSpecialData() {
    hat_mu_Y.fill(0);
    tilde_mu_X_prime.fill(0);
    lnDetV.fill(0);
    p.fill(0);
    Q.fill(0);
  }

  virtual void prepareBranch(uint i) {
    tTransf[i] = meta.t[i];
  }

  inline void pruneBranch(uint i) {
    if(i < meta.N) {
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

  inline void addToParent(uint i, uint iParent) {
    hat_mu_Y.row(iParent) += p[i]*hat_mu_Y.row(i);
    tilde_mu_X_prime.row(iParent) += p[i]*tilde_mu_X_prime.row(i);
    lnDetV[iParent] += lnDetV[i];
    p[iParent] += p[i];
    Q[iParent] += Q[i];
  }

  void do_pruning(int mode) {
    meta.do_pruning(this, mode);
  }

  const MetaType& get_meta() const {
    return meta;
  }
};
};
#endif // ParallelPruning_ThreePointV_lnDetV_Q_H_
