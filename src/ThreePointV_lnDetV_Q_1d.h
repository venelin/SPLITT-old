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

#ifndef ParallelPruning_ThreePointV_lnDetV_Q_1D_H_
#define ParallelPruning_ThreePointV_lnDetV_Q_1D_H_

#include "./ParallelPruningAlgorithm.h"
#include <iostream>

namespace ppa {
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
class ThreePointV_lnDetV_Q_1d {
protected:
  const ParallelPruningTree pptree;
  const ParallelPruningAlgorithm<ThreePointV_lnDetV_Q_1d> ppalgorithm;
  void init() {
    this->tTransf = vec(pptree.M - 1);
    this->lnDetV = vec(pptree.M, 0);
    this->p = vec(pptree.M, 0);
    this->Q = vec(pptree.M, 0);
  }
public:
  // define fields as public in order to access them easily from R.
  vec X, Y;
  vec tTransf;
  vec hat_mu_Y, tilde_mu_X_prime;
  vec lnDetV, p, Q;

  ThreePointV_lnDetV_Q_1d(Tree const& tree):
    pptree(tree), ppalgorithm(this->pptree, *this) {
    init();
  };

  void set_X_and_Y(vec const& X, vec const& Y) {
    if(X.size() != pptree.N || Y.size() != pptree.N) {
      Rcpp::stop("The matrices X and Y must have the same number of rows as V.");
    } else {
      this->X = X; this->Y = Y;

      this->hat_mu_Y = vec(pptree.M, 0);
      this->tilde_mu_X_prime = vec(pptree.M, 0);
    }
  }

  uint get_N() const {
    return pptree.get_N();
  }

  double get_Q() const {
    return this->Q[pptree.M-1];
  }

  double get_lnDetV() const {
    return this->lnDetV[pptree.M-1];
  }

  inline void initSpecialData() {
    std::fill(hat_mu_Y.begin(), hat_mu_Y.end(), 0);
    std::fill(tilde_mu_X_prime.begin(), tilde_mu_X_prime.end(), 0);
    std::fill(lnDetV.begin(), lnDetV.end(), 0);
    std::fill(p.begin(), p.end(), 0);
    std::fill(Q.begin(), Q.end(), 0);
  }

  virtual void prepareBranch(uint i) {
    tTransf[i] = pptree.t[i];
  }

  inline void pruneBranch(uint i) {
    if(i < pptree.N) {
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

  inline void addToParent(uint i, uint iParent) {
    hat_mu_Y[iParent] += p[i]*hat_mu_Y[i];
    tilde_mu_X_prime[iParent] += p[i]*tilde_mu_X_prime[i];
    lnDetV[iParent] += lnDetV[i];
    p[iParent] += p[i];
    Q[iParent] += Q[i];
  }
  void do_pruning(int mode) {
    ppalgorithm.do_pruning(mode);
  }
};
};

#endif // ParallelPruning_ThreePointV_lnDetV_Q_1D_H_
