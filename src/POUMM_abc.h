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
#ifndef POUMM_ABC_H_
#define POUMM_ABC_H_

#include "./ParallelPruningAlgorithm.h"

namespace ppa {
template<class Node>
class POUMM_abc {
  typedef ParallelPruningTree<Node, double> ParallelPruningTree;
  const ParallelPruningTree pptree;
  const ParallelPruningAlgorithm<ParallelPruningTree, POUMM_abc> ppalgorithm;
public:
  double alpha, theta, sigma, sigmae, sigmae2, sigma2, logsigma;
  vec z, se, a, b, c, sum_se2_sigmae2, talpha, etalpha,
  e2talpha, fe2talpha, gutalphasigma2;

  POUMM_abc(std::vector<Node> const& brStarts, std::vector<Node> const& brEnds, vec const& t,
            std::vector<Node> const& keys, vec const& z, vec const& se):
    pptree(brStarts, brEnds, t), ppalgorithm(this->pptree, *this), z(z), se(se) {

    if(z.size() != pptree.N || se.size() != pptree.N) {
      throw std::invalid_argument("The vectors z and se must be the same length as the number of tips.");
    } else {

      uvec ordNodes = pptree.order_nodes(keys);
      this->z = at(z, ordNodes);
      this->se = at(se, ordNodes);
      this->a = vec(pptree.M);
      this->b = vec(pptree.M);
      this->c = vec(pptree.M);
      this->talpha = vec(pptree.M - 1);
      this->etalpha = vec(pptree.M - 1);
      this->e2talpha = vec(pptree.M - 1);
      this->fe2talpha = vec(pptree.M - 1);
      this->gutalphasigma2 = vec(pptree.M - 1);
    }
  };

  vec get_abc() const {
    vec res(3);
    res[0] = a[pptree.M - 1];
    res[1] = b[pptree.M - 1];
    res[2] = c[pptree.M - 1];
    return res;
  };

  void set_parameters(double alpha, double theta, double sigma, double sigmae) {
    this->alpha = alpha;
    this->theta = theta;
    this->sigma = sigma;
    this->sigmae = sigmae;
    this->sigmae2 = sigmae*sigmae;

    this->sigma2 = sigma*sigma;
    this->logsigma = log(sigma);

    this->sum_se2_sigmae2 = se;
    for(size_t i = 0; i < se.size(); ++i) {
      sum_se2_sigmae2[i] = sigmae2 + se[i]*se[i];
    }
  };

  inline void initSpecialData() {
    std::fill(a.begin(), a.end(), 0);
    std::fill(b.begin(), b.end(), 0);
    std::fill(c.begin(), c.end(), 0);
  }

  inline void prepareBranch(uint i) {
    //std::cout<<"prepareBranch("<<i<<")"<<std::endl;
    if(alpha != 0) {
      talpha[i] = pptree.t[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = alpha / (1 - e2talpha[i]);
    } else {
      talpha[i] = pptree.t[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = -0.5 / pptree.t[i];
    }
  };

  inline void pruneBranch(uint i) {
    //std::cout<<"pruneBranch("<<i<<")"<<std::endl;
    if(i < pptree.N) {
      // branch leading to a tip
      gutalphasigma2[i] = e2talpha[i] +
        ((-0.5 / sum_se2_sigmae2[i]) * sigma2) / fe2talpha[i];
      double z1 = z[i] - theta;

      // integration over g1 including e1 = z1 - g1
      c[i] = -0.5 * log(gutalphasigma2[i]) -
        0.25 * sigma2 * z1*z1 / (sum_se2_sigmae2[i]*sum_se2_sigmae2[i]) /
          (fe2talpha[i] - alpha + (-0.5 / sum_se2_sigmae2[i]) * sigma2) +
            talpha[i] + (-0.5 * (M_LN_2PI  + z1*z1 / sum_se2_sigmae2[i]) -
            log(sqrt(sum_se2_sigmae2[i])));
      b[i] = (etalpha[i] * (z1 / sum_se2_sigmae2[i])) / gutalphasigma2[i];
      a[i] = (-0.5 / sum_se2_sigmae2[i]) / gutalphasigma2[i];
    } else {
      gutalphasigma2[i] = e2talpha[i] + (a[i] * sigma2) / fe2talpha[i];
      c[i] = -0.5 * log(gutalphasigma2[i]) - 0.25 * sigma2 * b[i] * b[i] /
        (fe2talpha[i] - alpha + a[i] * sigma2) + talpha[i] + c[i];
      b[i] = (etalpha[i] * b[i]) / gutalphasigma2[i];
      a[i] /= gutalphasigma2[i];
    }
  };

  inline void addToParent(uint i, uint iParent) {
    //std::cout<<"addToParent("<<i<<", "<<iParent<<")"<<std::endl;
    a[iParent] += a[i];
    b[iParent] += b[i];
    c[iParent] += c[i];
  }

  void do_pruning(int mode) {
    ppalgorithm.do_pruning(mode);
  }

};
};
#endif //POUMM_ABC_H_
