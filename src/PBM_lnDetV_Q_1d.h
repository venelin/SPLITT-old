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
#ifndef PBM_lnDetV_1D_H_
#define PBM_lnDetV_1D_H_

#include "ThreePointV_lnDetV_Q_1d.h"

namespace ppa {
class PBM_lnDetV_Q_1d: public ThreePointV_lnDetV_Q_1d {
public:
  // univariate trait vector
  vec z;
  double g0, sigma;
  PBM_lnDetV_Q_1d(Tree const& tree, vec const& z): ThreePointV_lnDetV_Q_1d(tree), z(z) {
    if(z.size() != meta.N ) {
      throw std::invalid_argument("The trait vector must have N elements (N is the number of tips).");
    } else {
      vec z2 = z;
      for(int i = 0; i < meta.N; ++i) {
        this->z[i] = z2[meta.orderNodes[i]];
      }
      vec X(meta.N);
      set_X_and_Y(X, X);
    }
  }

  void set_parameters(double g0, double sigma) {
    this->g0 = g0;
    this->sigma = sigma;
  }

  void prepareBranch(uint i) {
    tTransf[i] = meta.t[i]*sigma*sigma;
    if(i < meta.N) {
      X[i] = Y[i] = z[i] - g0;
    }
  }

  void do_pruning(int mode) {
    meta.do_pruning(this, mode);
  }

  const MetaType& get_meta() const {
    return meta;
  }
};

};

#endif // PBM_lnDetV_1D_H_
