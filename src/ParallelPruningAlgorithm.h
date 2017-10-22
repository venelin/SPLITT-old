// Copyright 2017 Venelin Mitov
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

// An Open-MP based implementation of a parallel pruning algorithm (postorder
// tree traversal). To use the library, follow the steps below:
// 0. The best way to familiarize youself with the five steps below is to learn
// the example class POUMM_abc.
//
// 1. Use a C++11 enabled compiler (option -std=c++11).
// 2. Include this header file in your code.
// 3. Create a class that inherits from ParallelPruningMeta:
// class MyPruningAlgorithm: public ParallelPruningMeta {...};
// 4. In the definition of your class, e.g. MyPruningAlgorithm:
//    4.1 add relevant data fields;
//    4.2 add code for setting the parameters for one pruning traversal
//    4.3 implement the following virtual methods:
//       4.3.1 A constructor MyPruningAlgorithm() that calls the parent class
//        constructor ParallelPruningMeta(N, branches_0, branches_1, t).
//       4.3.2 void prepareBranch(uint i)
//       4.3.3 void pruneBranch(uint i)
//       4.3.4 void addToParent(uint i, uint iParent)
//    4.4 If you allocated dynamic resources or memory, you might need to implement
//      a destructor ~MyPruningAlgorithm.
// 5. In your user code, create instances of your class, set the data fields and
//  and call do_pruning() on various parameter sets. After calling do_pruning(),
//  the pruning result at index M-1 should correspond to the final result at the
//  root of the tree.


#ifndef ParallelPruning_ParallelPruningAlgorithm_H_
#define ParallelPruning_ParallelPruningAlgorithm_H_

#include <algorithm>
#include <vector>
#include <math.h>
#include <iostream>
#include <sstream>
#include <limits>
#include <numeric>
#include <chrono>

// TODO: 1.  don't use Rcpp and arma classes but rely only on STL; done
// TODO: 2. create a wrapper interface with Rcpp; done
// TODO: 3. support for other tree descriptions (don't rely only on the phylo).
// TODO: 4. use std::exception to check arguments - done.

const uint MIN_CHUNK_SIZE = 10;

#ifdef _OPENMP

// Need to decide wheter to use '#pragma omp for' or '#pragma omp for simd'
#if _OPENMP >= 201307  // OMP 4.0 or higher

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#else // #if _OPENMP >= 201307

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // _OPENMP >= 201307

#include <omp.h>

#else // #ifdef _OPENMP

// the preprocessor directives should simply be ignored at compile-time
#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // #ifdef _OPENMP

// all functions and classes defined in the namespace ppa
// (stays for parallel pruning algorithm)
namespace ppa{

using namespace std;

typedef unsigned int uint;
typedef std::vector<uint> uvec;
typedef std::vector<double> vec;
typedef std::vector<bool> bvec;

// define an NA constant;
const uint NA_UINT = std::numeric_limits<uint>::max();

template <class VectorClass>
inline std::vector<size_t> order(VectorClass const& v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template<class VectorValues, class VectorPositions>
inline VectorValues at(VectorValues const& v, VectorPositions const& positions) {
  VectorValues sub;
  sub.resize(positions.size());

  size_t sub_i = 0;
  for(auto pit = positions.begin(); pit != positions.end(); pit++,sub_i++){
    sub[sub_i] = v[*pit];
  }
  return sub;
}

template<class VectorValues>
inline VectorValues at(VectorValues const& v, bvec const& mask) {
  if(mask.size() != v.size()) {
    throw std::length_error("bool vector mask should have the same length as v.");
  }

  size_t res_size = 0;
  for(auto b : mask) if(b) ++res_size;

  VectorValues sub(res_size);

  size_t sub_i = 0;
  for(uint i = 0; i < v.size(); i++){
    if(mask[i]) sub[sub_i++] = v[i];
  }
  return sub;
}

template<class VectorClass>
inline VectorClass multiReplace(VectorClass const& x, VectorClass const& a, VectorClass const& b) {
  auto x2 = x;
  auto ind = order(x2);
  auto xInd = at(x2, ind);
  std::pair<typename VectorClass::iterator, typename VectorClass::iterator> bounds;
  for(size_t i = 0; i < a.size(); ++i) {
    bounds = std::equal_range(xInd.begin(), xInd.end(), a.at(i));
    if(bounds.first != bounds.second) {
      size_t first = bounds.first - xInd.begin();
      size_t last = bounds.second - xInd.begin() - 1;
      for(size_t j = first; j<=last; ++j) {
        x2[ind[j]] = b[i];
      }
    }
  }
  return x2;
}


// for each element of x return the index of its first occurence in table or
// NA_UINT if the element is not found in table or is equal to NA_UINT.
// It is assumed that x does not have duplicated elements or NA_UINT elements.
inline uvec match(uvec const& x, uvec const& table) {
  auto minmax_x = std::minmax_element(x.begin(), x.end());
  uvec index(*minmax_x.second - *minmax_x.first + 1, NA_UINT);
  for(uint i = 0; i < table.size(); ++i) {
    if(table[i] >= *minmax_x.first && table[i] <= *minmax_x.second &&
       index[table[i] - *minmax_x.first] == NA_UINT) {
      index[table[i] - *minmax_x.first] = i;
    }
  }

  uvec positions(x.size());
  for(uint i = 0; i < x.size(); ++i) {
    positions[i] = index[x[i] - *minmax_x.first];
  }
  return positions;
}

inline uvec seq(uint first, uint last) {
  uvec res(last-first+1);
  std::iota(res.begin(), res.end(), first);
  return res;
}

inline bvec is_na(uvec const& x) {
  bvec res(x.size(), false);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA_UINT) res[i] = true;
  }
  return res;
}

inline bvec not_is_na(uvec const& x) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA_UINT) res[i] = false;
  }
  return res;
}


struct Tree {
  uint N;
  uvec branches_0, branches_1;
  vec t;
  Tree(uint N, uvec const& branches_0, uvec const& branches_1, vec const& t):
    N(N), branches_0(branches_0), branches_1(branches_1), t(t) {
    // check the parameters
    if(branches_0.size() != branches_1.size()) {
      std::ostringstream oss;
      oss<<"branches_0 and branches_1 should be the same size, but were "
         <<branches_0.size()<<" and "<<branches_1.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    if(t.size() != branches_0.size()) {
      std::ostringstream oss;
      oss<<"The size of t should be the same size as branches_0 and branches_1 ("
         <<branches_0.size()<<"), but was "<<t.size()<<".";
      throw std::length_error(oss.str());
    }
    if(N < 2) {
      std::ostringstream oss;
      oss<<"Trees with less than two tips are currently not supported. Check N.";
      throw std::length_error(oss.str());
    }
    if(*(min_element(branches_1.begin(), branches_1.end())) != 0) {
      std::ostringstream oss;
      oss<<"Tip 0 is not the min element in branches_1. "
         <<"The tips of the tree should be numerated from 0 to N-1.";
      throw std::invalid_argument(oss.str());
    }
    if(*(min_element(branches_0.begin(), branches_0.end())) != N) {
      std::ostringstream oss;
      oss<<"Node N should correspond to the root and be the min element in branches_0.";
      throw std::invalid_argument(oss.str());
    }
  }
};

class ParallelPruningTree {
  // call after set_tree()
  void createPruningOrder() {
    uvec branchEnds = branches_1;
    // insert a fictive branch leading to the root of the tree.
    branchEnds.push_back(N);

    uvec endingAt = match(seq(0, M - 1), branchEnds);


    uvec nonPrunedChildren(M);
    uvec ee1 = branches_0;

    auto start_ = std::chrono::steady_clock::now();
    while(ee1.size() > 0) {

      uvec matchp = match(seq(N, M - 1), ee1);

      matchp = at(matchp, not_is_na(matchp));

      for(auto m : matchp) nonPrunedChildren[ee1[m]]++;
      for(auto m : matchp) ee1[m] = NA_UINT;
      ee1 = at(ee1, not_is_na(ee1));

    }
    auto end_ = std::chrono::steady_clock::now();
    std::cout << "Duration: while(ee1...):" <<
      std::chrono::duration <double, milli> (end_ - start_).count() <<
        " ms" << std::endl;

    uvec tipsVector;
    uvec tipsVectorIndex(1, 0);

    uvec branchVector;
    uvec branchVectorIndex(1, 0);

    // start by pruning the tips
    uvec tips = seq(0, N-1);

    start_ = std::chrono::steady_clock::now();
    while(tips[0] != N) { // while the root has not become a tip itself

      tipsVectorIndex.push_back(
        tipsVectorIndex[tipsVectorIndex.size() - 1] + tips.size());


      // add the tips to be pruned to the tipsVector
      tipsVector.insert(tipsVector.end(), tips.begin(), tips.end());

      // indices in the branches-matrix of the to be pruned branches (pointing to tips)
      uvec branchesToTips = at(endingAt, tips);
      // empty the tips vector so it can be filled in with new tips to be pruned
      tips.clear();
      int nBranchesDone = 0;
      uvec remainingParents;

      while(nBranchesDone != branchesToTips.size() ) {
        // vectorized update of the parent nodes.
        // resolving order of sibling branches being pruned at the same time
        if(nBranchesDone == 0) {
          remainingParents = branches_0;
          remainingParents = at(remainingParents, branchesToTips);
        } else {
          remainingParents = branches_0;
          remainingParents = at(remainingParents,
                                at(branchesToTips, not_is_na(branchesToTips)));
        }


        //remainingParents = unique(remainingParents);
        std::sort(remainingParents.begin(), remainingParents.end());
        remainingParents = uvec(remainingParents.begin(),
                                std::unique(remainingParents.begin(),
                                            remainingParents.end()));


        uvec branchStarts(branchesToTips.size());

        for(int iett = 0; iett < branchesToTips.size(); iett++) {
          if(branchesToTips[iett] != NA_UINT) {
            branchStarts[iett] = branches_0[branchesToTips[iett]];
          } else {
            branchStarts[iett] = NA_UINT;
          }
        }

        // sib- branches that are first in branches[branchesToTips, 0] get served
        // first
        uvec branchesNext = match(remainingParents, branchStarts);
        std::sort(branchesNext.begin(), branchesNext.end());


        branchVector.insert(branchVector.end(),
                            branchesNext.begin(), branchesNext.end());

        // attach the index of the current last element of branchVector
        // to branchVectorIndex
        branchVectorIndex.push_back(
          branchVectorIndex[branchVectorIndex.size() - 1] + branchesNext.size());

        // For the parent nodes, decrement the amount of non-visited
        // children
        for(int u : branchesNext) {
          nonPrunedChildren[branches_0[branchesToTips[u]]]--;
          if(nonPrunedChildren[branches_0[branchesToTips[u]]] == 0) {
            tips.push_back(branches_0[branchesToTips[u]]);
          }
        }

        for(int u : branchesNext) branchesToTips[u] = NA_UINT;

        nBranchesDone += branchesNext.size();
      }
    }
    end_ = std::chrono::steady_clock::now();
    std::cout << "Duration: while(tips[0]!=N):" <<
      std::chrono::duration <double, milli> (end_ - start_).count() <<
        " ms" << std::endl;

    this->nLevels = tipsVectorIndex.size() - 1;
    this->tipsVector = tipsVector;
    this->tipsVectorIndex = tipsVectorIndex;
    this->branchVector = branchVector;
    this->branchVectorIndex = branchVectorIndex;
    this->orderNodes = seq(0, M - 1);

    start_ = std::chrono::steady_clock::now();

    reorderBranches();

    end_ = std::chrono::steady_clock::now();
    std::cout << "Duration: reorderBranches:" <<
      std::chrono::duration <double, milli> (end_ - start_).count() <<
        " ms" << std::endl;
  }

  void reorderBranches() {

    // the node number at which each branch ends
    uvec branchEnds(branches_1.size());
    copy(branches_1.begin(), branches_1.end(), branchEnds.begin());

    // add a branch ending at the root of the tree
    branchEnds.push_back(N);

    // the index of the branch ending at i; for the root node (N), endingAt[N]=(M-1)
    uvec endingAt(M);
    for(int i = 0; i < branchEnds.size(); ++i)
      endingAt[branchEnds[i]] = i;

    uvec parents = branches_0;
    replace(parents.begin(), parents.end(), N, 2*M - 1);

    (this->parentNode) = parents;
    // duplicate t and orderNodes since it will be shuffled during the reordering
    vec tOriginalOrder = (this->t);
    uvec orderNodesOriginal = (this->orderNodes);

    // branches pointing to tips
    uvec branchesToTips = at(
      endingAt, uvec(tipsVector.begin() + tipsVectorIndex[0],
                     tipsVector.begin() + tipsVectorIndex[1]));

    uint jBVI = 0;

    uint nBranchesDone = 0;
    auto start_ = std::chrono::steady_clock::now();
    while(nBranchesDone != branchesToTips.size()) {
      uvec branchesNext(branchVector.begin() + branchVectorIndex[jBVI],
                        branchVector.begin() + branchVectorIndex[jBVI + 1]);


      uvec branchEnds = at(branches_1, at(branchesToTips, branchesNext));

      uvec branchEndsNew(branchVectorIndex[jBVI + 1] - branchVectorIndex[jBVI]);
      std::iota(branchEndsNew.begin(),
                branchEndsNew.end(),
                branchVectorIndex[jBVI] + M);

      parents = multiReplace<uvec>(parents, branchEnds, branchEndsNew);

      auto parentNodeNew = at(parents, at(branchesToTips, branchesNext));
      std::copy(parentNodeNew.begin(), parentNodeNew.end(),
                (this->parentNode).begin()+branchVectorIndex[jBVI]);

      (this->parentNode) = multiReplace<uvec>(
        (this->parentNode), branchEnds, branchEndsNew);

      auto tNew = at(tOriginalOrder, at(branchesToTips, branchesNext));
      std::copy(tNew.begin(), tNew.end(),
                (this->t).begin() + branchVectorIndex[jBVI]);

      auto orderNodesNew = at(orderNodesOriginal, branchEnds);
      std::copy(orderNodesNew.begin(), orderNodesNew.end(),
                (this->orderNodes).begin() + branchVectorIndex[jBVI]);

      ++jBVI;
      nBranchesDone += branchesNext.size();
    }
    auto end_ = std::chrono::steady_clock::now();
    std::cout << "  Duration: while(nBranchesDone != branchesToTips.size()):" <<
      std::chrono::duration <double, milli> (end_ - start_).count() <<
        " ms" << std::endl;

    // branches pointing to internal nodes that have become tips
    for(int i = 1; i < nLevels; ++i) {
      branchesToTips = at(
        endingAt, uvec(tipsVector.begin() + tipsVectorIndex[i],
                       tipsVector.begin() + tipsVectorIndex[i + 1]));



      uint nBranchesDone = 0;
      while(nBranchesDone != branchesToTips.size()) {
        uvec branchesNext(branchVector.begin() + branchVectorIndex[jBVI],
                          branchVector.begin() + branchVectorIndex[jBVI + 1]);

        uvec branchEnds = at(branches_1, at(branchesToTips, branchesNext));

        uvec branchEndsNew(branchVectorIndex[jBVI + 1] - branchVectorIndex[jBVI]);
        std::iota(branchEndsNew.begin(),
                  branchEndsNew.end(),
                  branchVectorIndex[jBVI] + M);

        parents = multiReplace<uvec>(parents, branchEnds, branchEndsNew);

        auto parentNodeNew = at(parents, at(branchesToTips, branchesNext));
        std::copy(parentNodeNew.begin(), parentNodeNew.end(),
                  (this->parentNode).begin()+branchVectorIndex[jBVI]);

        (this->parentNode) = multiReplace<uvec>((this->parentNode), branchEnds, branchEndsNew);

        auto tNew = at(tOriginalOrder, at(branchesToTips, branchesNext));

        std::copy(tNew.begin(), tNew.end(),
                  (this->t).begin() + branchVectorIndex[jBVI]);

        auto orderNodesNew = at(orderNodesOriginal, branchEnds);

        std::copy(orderNodesNew.begin(), orderNodesNew.end(),
                  (this->orderNodes).begin() + branchVectorIndex[jBVI]);

        ++jBVI;
        nBranchesDone += branchesNext.size();
      }
    }

    for(int i = 0; i<parentNode.size(); ++i)
      (this->parentNode)[i] -= M;
  }

  // private default constructor;
  ParallelPruningTree() {};
public:
  uint nLevels = NA_UINT;
  // number of all nodes
  uint M;
  // number of tips
  uint N;

  // branches matrix
  uvec branches_0, branches_1;

  vec t;

  uvec parentNode;
  uvec orderNodes;

  uvec tipsVector;
  uvec tipsVectorIndex;
  uvec branchVector;
  uvec branchVectorIndex;

public:
  ParallelPruningTree(Tree const& tree): N(tree.N),
    branches_0(tree.branches_0), branches_1(tree.branches_1),
    t(tree.t) {
    // number of all nodes including the root
    this->M = branches_0.size() + 1;

    createPruningOrder();
  }

  uint get_nLevels() const {
    return nLevels;
  }

  uint get_M() const {
    return M;
  }

  uint get_N() const {
    return N;
  }

  vec get_t() const {
    return (t);
  }

  uvec get_parentNode() const {
    return (parentNode);
  }

  uvec get_tipsVector() const {
    return tipsVector;
  }

  uvec get_tipsVectorIndex() const {
    return tipsVectorIndex;
  }

  uvec get_branchVector() const {
    return branchVector;
  }

  uvec get_branchVectorIndex() const {
    return branchVectorIndex;
  }

  // The order in which the nodes get processed during the pruning, including
  // internal and tip-nodes.
  uvec get_orderNodes() const {
    return uvec(orderNodes.begin(), orderNodes.end() - 1);
  }

  // The order in which the branches get processed during the pruning.
  uvec get_orderBranches() const {
    uvec orderBranchEnds(orderNodes.begin(), orderNodes.end() - 1);
    return match(orderBranchEnds, branches_1);
  }

  // A root-to-node distance vector in the order of pruning processing
  vec get_nodeHeights() const {
    vec h(M, 0);
    for(int i = M - 2; i >= 0; i--) {
      h[i] = h[parentNode[i]] + t[i];
    }
    return h;
  }

};

//  a skeleton for a ParallelPruningSpec class
// Curiously recurring template pattern -  static polymorphism
template<class PruningSpec>
class ParallelPruningAlgorithm {
  uint nThreads;
protected:
  const ParallelPruningTree* p_tree;
  PruningSpec* p_spec;
public:
  ParallelPruningAlgorithm(ParallelPruningTree* p_tree, PruningSpec* p_spec):
  p_tree(p_tree), p_spec(p_spec) {
#ifdef _OPENMP
#pragma omp parallel
{
  uint tid = omp_get_thread_num();
  // only master thread does this
  if(tid == 0) {
    this->nThreads = omp_get_num_threads();
  }
}
#else
    this->nThreads = 1;
#endif // #ifdef _OPENMP
  }

  // void initSpecialData() {}
  //
  // void prepareBranch(uint i) {
  //   std::cout<<"preapareBranch("<<i<<")"<<std::endl;
  // };
  //
  // void pruneBranch(uint i) {
  //   std::cout<<"pruneBranch("<<i<<")"<<std::endl;
  // };
  //
  // void addToParent(uint i, uint iParent) {
  //   std::cout<<"addToParent("<<i<<","<<iParent<<")"<<std::endl;
  // };


  uint get_nThreads() const {
    return nThreads;
  }

  void do_pruning(int mode) const {
    switch(mode) {
    case 1: do_pruning_serial(); break;
    case 2: do_pruning_parallel(); break;
    default: do_pruning_hybrid();
    }
  }

  void do_pruning_serial() const {
    p_spec->initSpecialData();

    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < p_tree->M - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      p_spec->prepareBranch(i);
    }
    uint jBVI = 0;
    for(int j = 0; j < p_tree->nLevels; j++) {
      const uint bFirst = p_tree->tipsVectorIndex[j];
      const uint bLast = p_tree->tipsVectorIndex[j + 1] - 1;

      _PRAGMA_OMP_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        p_spec->pruneBranch(i);
      }
      uint nBranchesDone = 0;
      while(nBranchesDone != bLast - bFirst + 1) {
        const uint unFirst = p_tree->branchVectorIndex[jBVI];
        const uint unLast = p_tree->branchVectorIndex[jBVI + 1] - 1;

        _PRAGMA_OMP_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          p_spec->addToParent(i, p_tree->parentNode[i]);
        }

        nBranchesDone +=  p_tree->branchVectorIndex[jBVI + 1] - p_tree->branchVectorIndex[jBVI];
        ++jBVI;
      }
    }
  }

  void do_pruning_parallel() const {
    p_spec->initSpecialData();

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < p_tree->M - 1; i++) {
    // calculate and store cache information for branch i
    // This is the ideal place to do transformation of branch lengths, etc.
    p_spec->prepareBranch(i);
  }

  uint jBVI = 0;

  for(int j = 0; j < p_tree->nLevels; j++) {
    const uint bFirst = p_tree->tipsVectorIndex[j];
    const uint bLast = p_tree->tipsVectorIndex[j + 1] - 1;

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        p_spec->pruneBranch(i);
      }

      uint nBranchesDone = 0;

    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = p_tree->branchVectorIndex[jBVI];
      const uint unLast = p_tree->branchVectorIndex[jBVI + 1] - 1;

      _PRAGMA_OMP_FOR_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          p_spec->addToParent(i, p_tree->parentNode[i]);
        }

        nBranchesDone +=  p_tree->branchVectorIndex[jBVI + 1] - p_tree->branchVectorIndex[jBVI];
      ++jBVI;
#pragma omp barrier
    }
  }
}
  }

  void do_pruning_hybrid() const {
    p_spec->initSpecialData();

    uint tid;
#pragma omp parallel private(tid)
{
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  if(nThreads > 1) {
    _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < p_tree->M - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      p_spec->prepareBranch(i);
    }
  } else {
    // only one (master) thread executes this
    for(uint i = 0; i < p_tree->M - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      p_spec->prepareBranch(i);
    }
  }

  uint jBVI = 0;

  for(int j = 0; j < p_tree->nLevels; j++) {
    const uint bFirst = p_tree->tipsVectorIndex[j];
    const uint bLast = p_tree->tipsVectorIndex[j + 1] - 1;

    if(bLast - bFirst + 1 > nThreads * MIN_CHUNK_SIZE) {
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        p_spec->pruneBranch(i);
      }
    } else if(tid == 0) {
      // only one (master) thread executes this
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        p_spec->pruneBranch(i);
      }
    }

    uint nBranchesDone = 0;

    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = p_tree->branchVectorIndex[jBVI];
      const uint unLast = p_tree->branchVectorIndex[jBVI + 1] - 1;

      if(unLast - unFirst + 1 > nThreads * MIN_CHUNK_SIZE) {
        _PRAGMA_OMP_FOR_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          p_spec->addToParent(i, p_tree->parentNode[i]);
        }
      } else if(tid == 0) {
        // only one (master) thread executes this
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          p_spec->addToParent(i, p_tree->parentNode[i]);
        }
      }

      nBranchesDone +=  p_tree->branchVectorIndex[jBVI + 1] - p_tree->branchVectorIndex[jBVI];
      ++jBVI;
#pragma omp barrier
    }
  }
}
  }
};
}
#endif // ParallelPruning_ParallelPruningAlgorithm_H_
