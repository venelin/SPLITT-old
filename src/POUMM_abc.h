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
  const ParallelPruningTree pptree_;
  ParallelPruningAlgorithm<ParallelPruningTree, POUMM_abc> ppalgorithm_;
public:
  double alpha, theta, sigmae2, sigma2;
  vec z, se, a, b, c;

  POUMM_abc(std::vector<Node> const& brStarts, std::vector<Node> const& brEnds, vec const& t,
            std::vector<Node> const& keys, vec const& z, vec const& se):
    pptree_(brStarts, brEnds, t),
    ppalgorithm_(this->pptree_, *this),
    z(z), se(se) {

    if(z.size() != pptree_.num_tips() || se.size() != pptree_.num_tips()) {
      throw std::invalid_argument("The vectors z and se must be the same length as the number of tips.");
    } else {

      uvec ordNodes = pptree_.OrderNodes(keys);
      this->z = At(z, ordNodes);
      this->se = At(se, ordNodes);
      this->a = vec(pptree_.num_nodes());
      this->b = vec(pptree_.num_nodes());
      this->c = vec(pptree_.num_nodes());
    }
  };

  vec get_abc() const {
    vec res(3);
    res[0] = a[pptree_.num_nodes() - 1];
    res[1] = b[pptree_.num_nodes() - 1];
    res[2] = c[pptree_.num_nodes() - 1];
    return res;
  };

  void set_parameters(double alpha, double theta, double sigma, double sigmae) {
    this->alpha = alpha;
    this->theta = theta;
    this->sigmae2 = sigmae*sigmae;
    this->sigma2 = sigma*sigma;

  }

  uint num_threads() const {
    return ppalgorithm_.num_threads();
  }

  uint min_size_chunk_visit() const {
    return ppalgorithm_.min_size_chunk_visit();
  }

  uint min_size_chunk_prune() const {
    return ppalgorithm_.min_size_chunk_prune();
  }

  std::vector<double>  durations_tuning() const {
    return ppalgorithm_.durations_tuning();
  }

  uint ModeAuto() const {
    return ppalgorithm_.ModeAuto();
  }

  bool IsTuning() const {
    return ppalgorithm_.IsTuning();
  }
  uint IndexMinSizeChunkVisit() const {
    return ppalgorithm_.IndexMinSizeChunkVisit();
  }

  uint IndexMinSizeChunkPrune() const {
    return ppalgorithm_.IndexMinSizeChunkPrune();
  }

  uint fastest_step_tuning() const {
    return ppalgorithm_.fastest_step_tuning();
  }

  inline void InitNode(uint i) {
    if(i < this->pptree_.num_tips()) {
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

    double t = pptree_.LengthOfBranch(i);
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

  inline void PruneNode(uint i) {
    uint iParent = this->pptree_.FindIdOfParent(i);
    a[iParent] += a[i];
    b[iParent] += b[i];
    c[iParent] += c[i];
  }

  inline void DoPruning(int mode) {
    ppalgorithm_.DoPruning(mode);
  }
};

// typedef ParallelPruningInstance<
//   uint, double, ParallelPruningTree<uint, double>, AbcPOUMM::Data, AbcPOUMM> ParallelPruningAbcPOUMM;
//
//
}
#endif //POUMM_ABC_H_
