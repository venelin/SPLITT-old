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

// An Open-MP based implementation of a parallel pruning algorithm (postorder
// tree traversal). To use the library, follow the steps below:
// 0. The best way to familiarize youself with the five steps below is to learn
// the example class POUMM_abc.
//
// 1. Use a C++11 enabled compiler (option -std=c++11).
// 2. Include this header file in your code.
// 3. Create a class that inherits from ParallelPruningAlgorithm:
// class MyPruningAlgorithm: public ParallelPruningAlgorithm {...};
// 4. In the definition of your class, e.g. MyPruningAlgorithm:
//    4.1 add relevant data fields;
//    4.2 add code for setting the parameters for one pruning traversal
//    4.3 implement the following virtual methods:
//       4.3.1 A constructor MyPruningAlgorithm() that calls the parent class
//        constructor ParallelPruningAlgorithm(N, branches_0, branches_1, t).
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

// TODO: 1.  don't use Rcpp and arma classes but rely only on STL;
// TODO: 2. create a wrapper interface with Rcpp;
// TODO: 3. support for other tree descriptions (don't rely only on the phylo).
// TODO: 4. use std::exception to check arguments

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <R_ext/Rdynload.h>
#include <chrono>
#include <sstream>
#include <limits>


//[[Rcpp::plugins("cpp11")]]
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

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// all functions and classes defined in the namespace ppa
// (stays for parallel pruning algorithm)
namespace ppa{

using namespace std;

typedef unsigned int uint;
typedef std::vector<uint> uvec;
typedef std::vector<double> vec;
typedef std::vector<bool> bvec;

// define an NA constant;
const uint NA_uint = std::numeric_limits<uint>::max();

template <class VectorClass>
vector<size_t> order(VectorClass const& v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template<class VectorValues, class VectorPositions>
VectorValues at(VectorValues const& v, VectorPositions const& positions) {
  cout<<"    at:"<<endl;
  cout<<"    v:";
  for(auto m: v) cout<<m<<" ";
  cout<<endl;
  cout<<"    positions:";
  for(auto m: positions) cout<<m<<" ";
  cout<<endl;

  VectorValues sub;
  sub.resize(positions.size());

  size_t sub_i = 0;
  for(auto pit = positions.begin(); pit != positions.end(); pit++,sub_i++){
    sub[sub_i] = v[*pit];
  }
  return sub;
}

template<class VectorValues>
VectorValues at(VectorValues const& v, bvec const& mask) {
  cout<<"    at:"<<endl;
  cout<<"    v:";
  for(auto m: v) cout<<m<<" ";
  cout<<endl;
  cout<<"    mask:";
  for(auto m: mask) cout<<m<<" ";
  cout<<endl;

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
VectorClass multiReplace(VectorClass const& x, VectorClass const& a, VectorClass const& b) {
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
// NA_uint if the element is not found in table or is equal to NA_uint.
// It is assumed that x does not have duplicated elements or NA_uint elements.
uvec match(uvec const& x, uvec const& table) {
  cout<<"   match: "<<endl;

  cout<<"   x:";
  for(auto m: x) cout<<m<<" ";
  cout<<endl;

  cout<<"   table:";
  for(auto m: table) cout<<m<<" ";
  cout<<endl;

  auto minmax_x = std::minmax_element(x.begin(), x.end());
  uvec index(*minmax_x.second - *minmax_x.first + 1, NA_uint);
  for(uint i = 0; i < table.size(); ++i) {
    if(table[i] >= *minmax_x.first && table[i] <= *minmax_x.second &&
       index[table[i] - *minmax_x.first] == NA_uint) {
      index[table[i] - *minmax_x.first] = i;
    }
  }


  cout<<"   index:";
  for(auto m: index) cout<<m<<" ";
  cout<<endl;

  uvec positions(x.size());
  for(uint i = 0; i < x.size(); ++i) {
    positions[i] = index[x[i] - *minmax_x.first];
  }
  cout<<"   positions:";
  for(auto m: positions) cout<<m<<" ";
  cout<<endl;
  return positions;
}

uvec seq(uint first, uint last) {
  uvec res(last-first+1);
  std::iota(res.begin(), res.end(), first);
  return res;
}

bvec is_na(uvec const& x) {
  bvec res(x.size(), false);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA_uint) res[i] = true;
  }
  return res;
}

bvec not_is_na(uvec const& x) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA_uint) res[i] = false;
  }
  return res;
}


class ParallelPruningAlgorithm2 {
protected:
  uint nLevels;
  // number of all nodes
  uint M;
  // number of tips
  uint N;

  // branches matrix
  //arma::umat branches;
  uvec branches_0, branches_1;

  vec t;

  uvec parentNode;
  uvec orderNodes;

  uvec tipsVector;
  uvec tipsVectorIndex;
  uvec branchVector;
  uvec branchVectorIndex;

private:
  void set_inputTree(
      uint N, uvec const& branches_0, uvec const& branches_1, vec const& t) {

    cout<<"point1"<<endl;

    // check the parameters
    if(branches_0.size() != branches_1.size()) {
      std::ostringstream oss;
      oss<<"branches_0 and branches_1 should be the same size, but were "
      <<branches_0.size()<<" and "<<branches_1.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    cout<<"point2"<<endl;
    if(t.size() != branches_0.size()) {
      std::ostringstream oss;
      oss<<"The size of t should be the same size as branches_0 and branches_1 ("
         <<branches_0.size()<<"), but was "<<t.size()<<".";
      throw std::length_error(oss.str());
    }
    cout<<"point3"<<endl;
    if(N < 2) {
      std::ostringstream oss;
      oss<<"Trees with less than two tips are currently not supported. Check N.";
      throw std::length_error(oss.str());
    }
    cout<<"point4"<<endl;
    if(*(min_element(branches_1.begin(), branches_1.end())) != 0) {
      std::ostringstream oss;
      oss<<"Tip 0 is not the min element in branches_1. "
         <<"The tips of the tree should be numerated from 0 to N-1.";
      throw std::invalid_argument(oss.str());
    }
    cout<<"point5"<<endl;
    if(*(min_element(branches_0.begin(), branches_0.end())) != N) {
      std::ostringstream oss;
      oss<<"Node N should correspond to the root and be the min element in branches_0.";
      throw std::invalid_argument(oss.str());
    }
    cout<<"point6"<<endl;
    // assign parameters to fields.
    this->N = N;
    this->branches_0 = branches_0;
    this->branches_1 = branches_1;
    this->t = t;

    // number of all nodes including the root
    this->M = branches_0.size() + 1;
  }


  // call after set_inputTree()
  void createPruningOrder() {
    uvec branchEnds = branches_1;
    // insert a fictive branch leading to the root of the tree.
    branchEnds.push_back(N);

    uvec endingAt = match(seq(0, M - 1), branchEnds);

    cout<<"point7"<<endl;

    uvec nonPrunedChildren(M);
    uvec ee1 = branches_0;

    cout<<"point8"<<endl;
    while(ee1.size() > 0) {
      cout<<"point9"<<endl;

      cout<<"ee1:";
      for(auto m: ee1) cout<<m<<" ";
      cout<<endl;

      uvec matchp = match(seq(N, M - 1), ee1);

      cout<<"matchp:";
      for(auto m: matchp) cout<<m<<" ";
cout<<endl;
      matchp = at(matchp, not_is_na(matchp));



      for(auto m : matchp) nonPrunedChildren[ee1[m]]++;
      cout<<"point10"<<endl;
      for(auto m : matchp) ee1[m] = NA_uint;
      cout<<"point11"<<endl;
      ee1 = at(ee1, not_is_na(ee1));

      cout<<"point12"<<endl;
    }



    uvec tipsVector;
    uvec tipsVectorIndex(1, 0);

    uvec branchVector;
    uvec branchVectorIndex(1, 0);

    // start by pruning the tips
    uvec tips = seq(0, N-1);

    cout<<"point1"<<endl;
    while(tips[0] != N) { // while the root has not become a tip itself

      tipsVectorIndex.push_back(
        tipsVectorIndex[tipsVectorIndex.size() - 1] + tips.size());

      cout<<"point2"<<endl;

      // add the tips to be pruned to the tipsVector
      tipsVector.insert(tipsVector.end(), tips.begin(), tips.end());

      cout<<"point3"<<endl;
      // indices in the branches-matrix of the to be pruned branches (pointing to tips)
      uvec branchesToTips = at(endingAt, tips);
      cout<<"point4"<<endl;
      // empty the tips vector so it can be filled in with new tips to be pruned
      tips.clear();
      cout<<"point5"<<endl;
      int nBranchesDone = 0;
      uvec remainingParents;

      while(nBranchesDone != branchesToTips.size() ) {
        // vectorized update of the parent nodes.
        // resolving order of sibling branches being pruned at the same time
        if(nBranchesDone == 0) {
          remainingParents = branches_0;
          remainingParents = at(remainingParents, branchesToTips);
        } else {
          cout<<"point6.0"<<endl;
          remainingParents = branches_0;
          remainingParents = at(remainingParents,
                                at(branchesToTips, not_is_na(branchesToTips)));
        }

        cout<<"point6"<<endl;

        //remainingParents = unique(remainingParents);
        std::sort(remainingParents.begin(), remainingParents.end());
        remainingParents = uvec(remainingParents.begin(),
                                std::unique(remainingParents.begin(),
                                             remainingParents.end()));


        uvec branchStarts(branchesToTips.size());

        for(int iett = 0; iett < branchesToTips.size(); iett++) {
          if(branchesToTips[iett] != NA_uint) {
            branchStarts[iett] = branches_0[branchesToTips[iett]];
          } else {
            branchStarts[iett] = NA_uint;
          }
        }

        cout<<"point7"<<endl;
        // sib- branches that are first in branches[branchesToTips, 0] get served
        // first
        uvec branchesNext = match(remainingParents, branchStarts);
        std::sort(branchesNext.begin(), branchesNext.end());

        cout<<"point8"<<endl;

        branchVector.insert(branchVector.end(),
                            branchesNext.begin(), branchesNext.end());

        cout<<"point9"<<endl;
        // attach the index of the current last element of branchVector
        // to branchVectorIndex
        branchVectorIndex.push_back(
          branchVectorIndex[branchVectorIndex.size() - 1] + branchesNext.size());

        cout<<"point10"<<endl;
        // For the parent nodes, decrement the amount of non-visited
        // children
        for(int u : branchesNext) {
          nonPrunedChildren[branches_0[branchesToTips[u]]]--;
          if(nonPrunedChildren[branches_0[branchesToTips[u]]] == 0) {
            tips.push_back(branches_0[branchesToTips[u]]);
          }
        }

        cout<<"point11"<<endl;
        for(int u : branchesNext) branchesToTips[u] = NA_uint;

        cout<<"point12"<<endl;
        nBranchesDone += branchesNext.size();
      }
    }


    this->nLevels = tipsVectorIndex.size() - 1;
    this->tipsVector = tipsVector;
    this->tipsVectorIndex = tipsVectorIndex;
    this->branchVector = branchVector;
    this->branchVectorIndex = branchVectorIndex;
    this->orderNodes = seq(0, M - 1);
    cout<<"point13"<<endl;
    reorderBranches();
  }

  void createPruningOrder(Rcpp::List const& tree) {
    using namespace Rcpp;
    if(!tree.inherits("phylo")) {
      Rcpp::stop("Input must be a phylo object.");
    }

    // 0-based indices
    Rcpp::IntegerMatrix branches = Rcpp::as<IntegerMatrix>(tree["edge"]) - 1;

    IntegerVector br_0 = branches(_, 0);
    IntegerVector br_1 = branches(_, 1);

    NumericVector branchLength = tree["edge.length"];

    set_inputTree(N = Rcpp::as<CharacterVector>(tree["tip.label"]).size(),
                  branches_0 = uvec(br_0.begin(), br_0.end()),
                  branches_1 = uvec(br_1.begin(), br_1.end()),
                  t = Rcpp::as<vec>(branchLength));

    createPruningOrder();

    // // number of tips
    // this->N = Rcpp::as<CharacterVector>(tree["tip.label"]).size();
    // // number of all nodes
    // this->M = unique(Rcpp::as<IntegerVector>(branches)).size();
    //
    // IntegerVector branchEnds = branches(_, 1);
    // branchEnds.push_back(N);
    //
    // IntegerVector endingAt = Rcpp::match(Rcpp::seq(0, M - 1), branchEnds) - 1;
    //
    // IntegerVector nonPrunedChildren(M);
    // IntegerVector ee1 = branches(_, 0);
    //
    // while(ee1.size() > 0) {
    //   IntegerVector matchp = Rcpp::match(Rcpp::seq(N, M - 1), ee1) - 1;
    //   matchp = matchp[!Rcpp::is_na(matchp)];
    //
    //   for(int m : matchp) {
    //     nonPrunedChildren[ee1[m]]++;
    //   }
    //
    //   ee1[matchp] = -1;
    //   ee1 = ee1[ee1 != (-1)];
    // }
    //
    // uvec tipsVector;
    // uvec tipsVectorIndex(1, 0);
    //
    // uvec branchVector;
    // uvec branchVectorIndex(1, 0);
    //
    // // start by pruning the tips
    // uvec tips = seq(0, N-1);
    //
    // while(tips[0] != N) { // while the root has not become a tip itself
    //
    //   tipsVectorIndex.push_back(
    //     tipsVectorIndex[tipsVectorIndex.size() - 1] + tips.size());
    //
    //   // add the tips to be pruned to the tipsVector
    //   tipsVector.insert(tipsVector.end(), tips.begin(), tips.end());
    //
    //   // indices in the branches-matrix of the to be pruned branches (pointing to tips)
    //   IntegerVector branchesToTips = endingAt[IntegerVector(tips.begin(), tips.end())];
    //
    //   // empty the tips vector so it can be filled in with new tips to be pruned
    //   tips.clear();
    //
    //   int nBranchesDone = 0;
    //   IntegerVector remainingParents;
    //
    //   while(nBranchesDone != branchesToTips.size() ) {
    //     // vectorized update of the parent nodes.
    //     // resolving order of sibling branches being pruned at the same time
    //     if(nBranchesDone == 0) {
    //       remainingParents = branches(_, 0);
    //       remainingParents = remainingParents[branchesToTips];
    //     } else {
    //       remainingParents = branches(_, 0);
    //       remainingParents = remainingParents[
    //       Rcpp::as<IntegerVector>(branchesToTips[! is_na(branchesToTips)])];
    //     }
    //
    //     remainingParents = unique(remainingParents);
    //
    //     IntegerVector branchStarts(branchesToTips.size());
    //     for(int iett = 0; iett < branchesToTips.size(); iett++) {
    //       if(branchesToTips(iett) != NA_INTEGER) {
    //         branchStarts(iett) = branches(branchesToTips(iett), 0);
    //       } else {
    //         branchStarts(iett) = NA_INTEGER;
    //       }
    //     }
    //
    //     // sib- branches that are first in branches[branchesToTips, 0] get served
    //     // first
    //     IntegerVector branchesNext = Rcpp::match(remainingParents, branchStarts);
    //     branchesNext = branchesNext - 1;
    //     branchesNext = branchesNext.sort();
    //
    //     branchVector.insert(branchVector.end(), branchesNext.begin(), branchesNext.end());
    //
    //     // attach the index of the current last element of branchVector
    //     // to branchVectorIndex
    //     branchVectorIndex.push_back(
    //       branchVectorIndex[branchVectorIndex.size() - 1] + branchesNext.size());
    //
    //     // For the parent nodes, decrement the amount of non-visited
    //     // children
    //     for(int u : branchesNext) {
    //       nonPrunedChildren[branches(branchesToTips[u], 0)]--;
    //       if(nonPrunedChildren[branches(branchesToTips[u], 0)] == 0) {
    //         tips.push_back(branches(branchesToTips[u], 0));
    //       }
    //     }
    //
    //     for(int u : branchesNext) branchesToTips[u] = NA_INTEGER;
    //
    //     nBranchesDone += branchesNext.size();
    //   }
    // }
    //
    //
    // this->nLevels = tipsVectorIndex.size() - 1;
    //
    //
    // this->tipsVector = tipsVector;
    //
    // this->tipsVectorIndex = tipsVectorIndex;
    //
    // this->branchVector = branchVector;
    //
    // this->branchVectorIndex = branchVectorIndex;
    //
    // IntegerVector br_0 = branches.column(0);
    // this->branches_0 = as<uvec>(br_0);
    //
    // IntegerVector br_1 = branches.column(1);
    // this->branches_1 = as<uvec>(br_1);
    //
    // NumericVector branchLength = tree["edge.length"];
    // this->t = Rcpp::as<vec>(branchLength);
    //
    // IntegerVector orderNodesOriginal = Rcpp::seq(0, M - 1);
    // this->orderNodes = Rcpp::as<uvec>(orderNodesOriginal);
//
    // reorderBranches();
  };

  void reorderBranches() {

    uvec branchEnds(branches_1.size());
    copy(branches_1.begin(), branches_1.end(), branchEnds.begin());
    branchEnds.push_back(N);

    uvec endingAt(M);
    for(int i = 0; i < branchEnds.size(); ++i)
      endingAt[branchEnds[i]] = i;

    cout<<"point14"<<endl;

    uvec parents = branches_0;
    replace(parents.begin(), parents.end(), N, 2*M - 1);

    (this->parentNode) = parents;
    // duplicate t and orderNodes since it will be shuffled during the reordering
    vec tOriginalOrder = (this->t);
    uvec orderNodesOriginal = (this->orderNodes);

    cout<<"point15"<<endl;

    // branches pointing to tips
    uvec branchesToTips = at(
      endingAt, uvec(tipsVector.begin() + tipsVectorIndex[0],
                     tipsVector.begin() + tipsVectorIndex[1]));

    cout<<"point16"<<endl;

    uint jBVI = 0;

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
  };

public:
  // Default constructor;
  ParallelPruningAlgorithm2(Rcpp::List const& tree) {
    this->createPruningOrder(tree);
  };

  uint get_nLevels() const {
    return nLevels;
  }

  uint get_M() const {
    return M;
  }

  uint get_N() const {
    return N;
  }

  arma::vec get_t() const {
    return (this->t);
  }

  arma::uvec get_parentNode() const {
    return (this->parentNode);
  }

  arma::uvec get_tipsVector() const {
    return tipsVector;
  }

  arma::uvec get_tipsVectorIndex() const {
    return tipsVectorIndex;
  }

  arma::uvec get_branchVector() const {
    return branchVector;
  }

  arma::uvec get_branchVectorIndex() const {
    return branchVectorIndex;
  }


  virtual void set_parameters(Rcpp::List const& par) {
    cout<<"setParameters()"<<endl;
  };

  virtual void prepareBranch(uint i) {
    cout<<"preapareCachForBranch("<<i<<")"<<endl;
  };

  virtual void pruneBranch(uint i) {
    cout<<"pruneBranch("<<i<<")"<<endl;
  };
  virtual void addToParent(uint i, uint iParent) {
    cout<<"addToParent("<<i<<","<<iParent<<")"<<endl;
  };

  void do_pruning(Rcpp::List const& par) {
    set_parameters(par);

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < M - 1; i++) {
    // calculate and store cache information for branch i
    prepareBranch(i);
  }
  uint jBVI = 0;

  for(int j = 0; j < nLevels; j++) {
    const uint bFirst = tipsVectorIndex[j];
    const uint bLast = tipsVectorIndex[j + 1] - 1;

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        pruneBranch(i);
      }

      uint nBranchesDone = 0;
    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = branchVectorIndex[jBVI];
      const uint unLast = branchVectorIndex[jBVI + 1] - 1;
      _PRAGMA_OMP_FOR_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          addToParent(i, parentNode[i]);
        }
        nBranchesDone +=  branchVectorIndex[jBVI + 1] - branchVectorIndex[jBVI];
      ++jBVI;
    }
  }
}
  };

  virtual ~ParallelPruningAlgorithm2() {};
};

RCPP_MODULE(ParallelPruningAlgorithm2) {
  Rcpp::class_<ParallelPruningAlgorithm2>( "ParallelPruningAlgorithm2" )
  .constructor<Rcpp::List const&>()
  .method( "do_pruning", &ParallelPruningAlgorithm2::do_pruning )
  .property("nLevels", &ParallelPruningAlgorithm2::get_nLevels )
  .property("M", &ParallelPruningAlgorithm2::get_M )
  .property("N", &ParallelPruningAlgorithm2::get_N )
  .property("t", &ParallelPruningAlgorithm2::get_t )
  .property("parentNode", &ParallelPruningAlgorithm2::get_parentNode )
  .property("tipsVector", &ParallelPruningAlgorithm2::get_tipsVector )
  .property("tipsVectorIndex", &ParallelPruningAlgorithm2::get_tipsVectorIndex )
  .property("branchVector", &ParallelPruningAlgorithm2::get_branchVector )
  .property("branchVectorIndex", &ParallelPruningAlgorithm2::get_branchVectorIndex )
  ;
}

class POUMM_abc: public ParallelPruningAlgorithm2 {
protected:
  arma::vec z;
  arma::vec se;

  bool has_se;

  arma::vec sum_se2_sigmae2;

  arma::mat abcMat;
  arma::vec a, b, c;

  uint count_abc_calls;


  double g0, alpha, theta, sigma, sigmae;
  double sigmae2, sigma2, logsigma;

  arma::vec log_se_total;
  arma::vec talpha;
  arma::vec etalpha;
  arma::vec e2talpha;
  arma::vec fe2talpha;
  arma::vec gutalphasigma2;

  arma::vec z1, z1z1;

public:
  POUMM_abc(Rcpp::List const& tree,
            Rcpp::NumericVector const& z,
            Rcpp::NumericVector const& se): ParallelPruningAlgorithm2(tree) {
    if(z.size() != N) {
      Rcpp::stop("The trait vector z must be the same length as the number of tips.");
    } else {
      arma::vec zArma = z;
      arma::vec seArma = se;
      this->z = zArma;
      this->se = seArma;

      for(int i = 0; i < N; ++i) {
        this->z(i) = zArma(orderNodes[i]);
        this->se(i) = seArma(orderNodes[i]);
      }

      this->a = arma::vec(M);
      this->b = arma::vec(M);
      this->c = arma::vec(M);
      this->talpha = arma::vec(M - 1);
      this->etalpha = arma::vec(M - 1);
      this->e2talpha = arma::vec(M - 1);
      this->fe2talpha = arma::vec(M - 1);
      this->gutalphasigma2 = arma::vec(M - 1);
      this->z1 = arma::vec(N);
      this->z1z1 = arma::vec(N);
    }
  };

  arma::vec get_abc() const {
    arma::vec res(3);
    res(0) = a(M - 1);
    res(1) = b(M - 1);
    res(2) = c(M - 1);
    return res;
  };

  arma::mat get_abcMat() const {
    arma::mat res(M, 3);
    res.col(0) = this->a;
    res.col(1) = this->b;
    res.col(2) = this->c;
    return res;
  };

  void set_parameters(Rcpp::List const& par) {
    this->g0 = par["g0"];
    this->alpha = par["alpha"];
    this->theta = par["theta"];
    this->sigma = par["sigma"];
    this->sigmae = par["sigmae"];

    this->sigmae2 = sigmae*sigmae;
    this->sum_se2_sigmae2 = sigmae2 + se % se;

    this->log_se_total = log(sqrt(sum_se2_sigmae2));

    this->sigma2 = sigma*sigma;
    this->logsigma = log(sigma);

    this->a.fill(0);
    this->b.fill(0);
    this->c.fill(0);
  };

  inline void prepareBranch(uint i) {
    if(alpha != 0) {
      talpha[i] = (this->t)[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = alpha / (1 - e2talpha[i]);
    } else {
      talpha[i] = (this->t)[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = -0.5 / (this->t)[i];
    }
  };

  inline void pruneBranch(uint i) {
    if(i < N) {
      // branch leading to a tip
      double z1, z1z1;
      gutalphasigma2[i] = e2talpha[i] + ((-0.5 / sum_se2_sigmae2[i]) * sigma2) / fe2talpha[i];
      z1 = z[i] - theta;
      z1z1 = z1 * z1;

      // integration over g1 including e1 = z1 - g1
      c[i] = -0.5 * log(gutalphasigma2[i]) -
        0.25 * sigma2 * z1z1 / (sum_se2_sigmae2[i]*sum_se2_sigmae2[i]) /
          (fe2talpha[i] - alpha + (-0.5 / sum_se2_sigmae2[i]) * sigma2) +
            talpha[i] + (-0.5 * (M_LN_2PI  + z1z1 / sum_se2_sigmae2[i]) - log_se_total[i]);
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
    a[iParent] += a[i];
    b[iParent] += b[i];
    c[iParent] += c[i];
  }
};

RCPP_MODULE(POUMM_abc) {
  Rcpp::class_<ParallelPruningAlgorithm2>( "ParallelPruningAlgorithm2" )
    .constructor<Rcpp::List const&>()
    .method( "do_pruning", &ParallelPruningAlgorithm2::do_pruning )
    .property("nLevels", &ParallelPruningAlgorithm2::get_nLevels )
    .property("M", &ParallelPruningAlgorithm2::get_M )
    .property("N", &ParallelPruningAlgorithm2::get_N )
    .property("t", &ParallelPruningAlgorithm2::get_t )
    .property("parentNode", &ParallelPruningAlgorithm2::get_parentNode )
    .property("tipsVector", &ParallelPruningAlgorithm2::get_tipsVector )
    .property("tipsVectorIndex", &ParallelPruningAlgorithm2::get_tipsVectorIndex )
    .property("branchVector", &ParallelPruningAlgorithm2::get_branchVector )
    .property("branchVectorIndex", &ParallelPruningAlgorithm2::get_branchVectorIndex )
  ;
  Rcpp::class_<POUMM_abc>( "POUMM_abc" )
    .derives<ParallelPruningAlgorithm2>("ParallelPruningAlgorithm2")
    .constructor<Rcpp::List const&, Rcpp::NumericVector const&, Rcpp::NumericVector const&>()
    .method( "abcMat", &POUMM_abc::get_abcMat )
    .method( "abc", &POUMM_abc::get_abc )
  ;
}
}

#endif ParallelPruning_ParallelPruningAlgorithm_H_
