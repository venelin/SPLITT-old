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
#include <map>

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
inline std::vector<uint> order(VectorClass const& v) {

  // initialize original index locations
  std::vector<uint> idx(v.size());
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

// A tree topology is a tree without branch-lengths. Still, for our purposes,
// it has named nodes. Permuting the nodes results in a different topology.
template<class Node, class BranchWeight>
class Tree {
protected:
  uint N;
  uint M;
  uvec branches_0, branches_1;
  std::map<Node, uint> mapNodeToId;
  std::vector<Node> mapIdToNode;
  std::vector<BranchWeight> t;

public:
  Tree(std::vector<Node> const& brStarts,
       std::vector<Node> const& brEnds,
       std::vector<BranchWeight> const& t): t(t) {
    if(brStarts.size() != brEnds.size()) {
      std::ostringstream oss;
      oss<<"brStarts and brEnds should be the same size, but were "
         <<brStarts.size()<<" and "<<brEnds.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    // There should be exactly M = number-of-branches+1 distinct nodes
    // This is because each branch can be mapped to its ending node. The +1
    // corresponds to the root node, to which no branch points.
    this->M = brStarts.size() + 1;

    // we distinguish three types of nodes:
    enum NodeType { ROOT, INTERNAL, TIP };

    // initially, we traverse the list of branches and order the nodes as they
    // appear during the traversal from 0 to M-1. This order is done by
    // incrementing nodeIdTemp.
    uint nodeIdTemp = 0;

    std::vector<NodeType> nodeTypes(M, ROOT);
    this->mapIdToNode = std::vector<Node>(M);
    uvec branchStartsTemp(brStarts.size(), NA_UINT);
    uvec branchEndsTemp(brStarts.size(), NA_UINT);
    uvec endingAt(M - 1, NA_UINT);

    for(uint i = 0; i < brStarts.size(); ++i) {
      if(brStarts[i] == brEnds[i]) {
        std::ostringstream oss;
        oss<<"Found a branch with the same start and end node ("<<
          brStarts[i]<<"). Not allowed. ";
        throw std::logic_error(oss.str());
      }

      auto it1 = mapNodeToId.find(brStarts[i]);
      if(it1 == mapNodeToId.end()) {
        // node encountered for the first time
        mapNodeToId[brStarts[i]] = nodeIdTemp; // insert brStarts[i] in the map
        mapIdToNode[nodeIdTemp] = brStarts[i];
        if(nodeTypes[nodeIdTemp] == TIP) {
          nodeTypes[nodeIdTemp] = INTERNAL;
        }
        branchStartsTemp[i] = nodeIdTemp;
        nodeIdTemp++;
      } else {
        // node encountered in a previous branch
        if(nodeTypes[it1->second] == TIP) {
          // the previous encounter of the node was as a branch-end
          nodeTypes[it1->second] = INTERNAL;
        } else {
          // do nothing
        }
        branchStartsTemp[i] = it1->second;
      }

      auto it2 = mapNodeToId.find(brEnds[i]);
      if(it2 == mapNodeToId.end()) {
        // node encountered for the first time
        mapNodeToId[brEnds[i]] = nodeIdTemp;
        mapIdToNode[nodeIdTemp] = brEnds[i];

        if(nodeTypes[nodeIdTemp] == ROOT) {
          // not known if the node has descendants, so we set its type to TIP.
          nodeTypes[nodeIdTemp] = TIP;
        }
        branchEndsTemp[i] = nodeIdTemp;
        endingAt[nodeIdTemp] = i;
        nodeIdTemp++;
      } else {
        // node has been previously encountered
        if(endingAt[it2->second] != NA_UINT) {
          std::ostringstream oss;
          oss<<"Found at least two branches ending at the same node ("<<
            it2->first<<"). Check for cycles or repeated branches. ";
          throw std::logic_error(oss.str());
        } else {
          if(nodeTypes[it2->second] == ROOT) {
            // the previous enounters of the node were as branch-start -> set
            // the node's type to INTERNAL, because we know for sure that it
            // has descendants.
            nodeTypes[it2->second] = INTERNAL;
          }
          branchEndsTemp[i] = it2->second;
          endingAt[it2->second] = i;
        }
      }
    }

    if(mapNodeToId.size() != M) {
      std::ostringstream oss;
      oss<<"The number of distinct nodes ("<<mapNodeToId.size()<<
        ") should equal the number-of-branches+1 ("<<M<<").";
      throw std::logic_error(oss.str());
    }

    auto countRoots = count(nodeTypes.begin(), nodeTypes.end(), ROOT);
    if(countRoots != 1) {
      std::ostringstream oss;
      oss<<"There should be exactly one ROOT node, but "<<countRoots<<
        " were found. Check for cycles or for multiple trees.";
      throw std::logic_error(oss.str());
    }

    auto countTips = count(nodeTypes.begin(), nodeTypes.end(), TIP);
    if(countTips == 0) {
      std::ostringstream oss;
      oss<<"There should be at least one TIP node, but none"<<
        " was found. Check for cycles.";
      throw std::logic_error(oss.str());
    }
    this->N = countTips;

    // assign new ids according to the following convention:
    // tips are numbered from 0 to N - 1;
    // internal nodes are numbered from N to M - 2;
    // root is numbered M - 1;
    std::vector<uint> nodeIds(M, NA_UINT);

    uint tipNo = 0, internalNo = N, rootNo = M - 1;
    for(uint i = 0; i < M; i++) {
      if(nodeTypes[i] == TIP) {
        nodeIds[i] = tipNo;
        tipNo ++;
      } else if(nodeTypes[i] == INTERNAL) {
        nodeIds[i] = internalNo;
        internalNo++;
      } else {
        // nodeTypes[i] == ROOT
        nodeIds[i] = M - 1;
      }
      mapNodeToId[mapIdToNode[i]] = nodeIds[i];
    }

    this->mapIdToNode = at(mapIdToNode, order(nodeIds));

    this->branches_0 = uvec(M - 1);
    this->branches_1 = uvec(M - 1);

    for(uint i = 0; i < M - 1; i++) {
      branches_0[i] = nodeIds[branchStartsTemp[i]];
      branches_1[i] = nodeIds[branchEndsTemp[i]];
    }
  }

  uint get_M() const {
    return M;
  }

  uint get_N() const {
    return N;
  }

  vec get_t() const {
    return t;
  }

  // return the node at position id;
  Node get_node(uint id) const {
    return mapIdToNode[id];
  }

  // returns the vector of nodes in order of their ids:
  // tips at positions from 0 to N-1;
  // internal nodes ath positions from N to M-2;
  // rot at position M-1;
  vector<Node> get_nodes() const {
    return mapIdToNode;
  }

  // get the internally stored id of node
  uint get_id(Node const& node) const {
    return mapNodeToId[node];
  }
};

template<class Node, class BranchWeight>
class ParallelPruningTree: public Tree<Node, BranchWeight> {

  void createPruningOrder() {
    uvec branchEnds = this->branches_1;
    // insert a fictive branch leading to the root of the tree.
    branchEnds.push_back(this->M - 1);

    uvec endingAt = match(seq(0, this->M - 1), branchEnds);


    uvec nonPrunedChildren(this->M);
    uvec ee1 = this->branches_0;

    while(ee1.size() > 0) {
      uvec matchp = match(seq(this->N, this->M - 1), ee1);
      matchp = at(matchp, not_is_na(matchp));
      for(auto m : matchp) nonPrunedChildren[ee1[m]]++;
      for(auto m : matchp) ee1[m] = NA_UINT;
      ee1 = at(ee1, not_is_na(ee1));

    }

    uvec tipsVector;
    uvec tipsVectorIndex(1, 0);

    uvec branchVector;
    uvec branchVectorIndex(1, 0);

    // start by pruning the tips
    uvec tips = seq(0, this->N-1);

    while(tips[0] != this->M-1) { // while the root has not become a tip itself

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
          remainingParents = this->branches_0;
          remainingParents = at(remainingParents, branchesToTips);
        } else {
          remainingParents = this->branches_0;
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
            branchStarts[iett] = this->branches_0[branchesToTips[iett]];
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
          nonPrunedChildren[this->branches_0[branchesToTips[u]]]--;
          if(nonPrunedChildren[this->branches_0[branchesToTips[u]]] == 0) {
            tips.push_back(this->branches_0[branchesToTips[u]]);
          }
        }

        for(int u : branchesNext) branchesToTips[u] = NA_UINT;

        nBranchesDone += branchesNext.size();
      }
    }

    this->nLevels = tipsVectorIndex.size() - 1;
    this->tipsVector = tipsVector;
    this->tipsVectorIndex = tipsVectorIndex;
    this->branchVector = branchVector;
    this->branchVectorIndex = branchVectorIndex;
    this->orderNodeIds = seq(0, this->M - 1);

    reorderBranches(branchEnds, endingAt);
  }

  void reorderBranches(uvec const& branchEnds, uvec const& endingAt) {
    // duplicate orderNodeIds since it will be shuffled during the reordering
    uvec orderNodesOriginal = (this->orderNodeIds);

    uint jBVI = 0;
    uint nBranchesDone = 0;

    // branches pointing to internal nodes that have become tips
    for(int i = 0; i < nLevels; ++i) {
      uvec branchesToTips = at(
        endingAt, uvec(tipsVector.begin() + tipsVectorIndex[i],
                       tipsVector.begin() + tipsVectorIndex[i + 1]));

      uint nBranchesDone = 0;
      while(nBranchesDone != branchesToTips.size()) {
        uvec branchesNext(branchVector.begin() + branchVectorIndex[jBVI],
                          branchVector.begin() + branchVectorIndex[jBVI + 1]);

        uvec branchEnds = at(this->branches_1, at(branchesToTips, branchesNext));

        auto orderNodesNew = at(orderNodesOriginal, branchEnds);
        std::copy(orderNodesNew.begin(), orderNodesNew.end(),
                  (this->orderNodeIds).begin() + branchVectorIndex[jBVI]);

        ++jBVI;
        nBranchesDone += branchesNext.size();
      }
    }

    auto orderBranches = get_orderBranches();
    this->t = at(this->t, orderBranches);

    uvec branchEndsReord = at(this->branches_1, orderBranches);
    branchEndsReord.push_back(this->M - 1);
    this->parentNode = match(at(this->branches_0, orderBranches), branchEndsReord);
  }

  // private default constructor;
  ParallelPruningTree() {};
public:
  // number of pruning steps
  uint nLevels;

  uvec parentNode;
  uvec orderNodeIds;

  uvec tipsVector;
  uvec tipsVectorIndex;
  uvec branchVector;
  uvec branchVectorIndex;

public:
  ParallelPruningTree(
    std::vector<Node> const& brStarts,
    std::vector<Node> const& brEnds,
    std::vector<BranchWeight> const& t): Tree<Node, BranchWeight>(brStarts, brEnds, t) {

    createPruningOrder();
  }

  uint get_nLevels() const {
    return nLevels;
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

  // The order in which the node ids get processed during the pruning, including
  // tips, internal and root nodes. Each element of the returned
  // vector is an internal node-id. See also get_orderNodes.
  uvec get_orderNodeIds() const {
    return uvec(orderNodeIds.begin(), orderNodeIds.end());
  }

  // The order in which the branches get processed during the pruning.
  uvec get_orderBranches() const {
    uvec orderBranchEnds(orderNodeIds.begin(), orderNodeIds.end() - 1);
    return match(orderBranchEnds, this->branches_1);
  }

  // A root-to-node distance vector in the order of pruning processing
  vec get_nodeHeights() const {
    vec h(this->M, 0);
    for(int i = this->M - 2; i >= 0; i--) {
      h[i] = h[parentNode[i]] + this->t[i];
    }
    return h;
  }

  // returns a vector of positions in nodes in the order of pruning.
  // See also get_orderNodeIds.
  uvec order_nodes(std::vector<Node> const& nodes) const {
    uvec ordIds(nodes.size());
    for(uint i = 0; i < nodes.size(); ++i) {
      auto it = this->mapNodeToId.find(nodes[i]);
      if(it == this->mapNodeToId.end()) {
        std::ostringstream oss;
        oss<<"At least one of the nodes is not present in the tree.";
        throw std::invalid_argument(oss.str());
      } else {
        ordIds[i] = orderNodeIds[it->second];
      }
    }
    return order(ordIds);
  }
};

//  a skeleton for a ParallelPruningSpec class
// Curiously recurring template pattern -  static polymorphism
template<class ParallelPruningTree, class PruningSpec>
class ParallelPruningAlgorithm {
  uint nThreads;
protected:
  const ParallelPruningTree& pptree;
  PruningSpec& spec;
public:
  ParallelPruningAlgorithm(const ParallelPruningTree& _pptree, PruningSpec& _spec):
  pptree(_pptree), spec(_spec) {
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
    spec.initSpecialData();

    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < pptree.M - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      spec.prepareBranch(i);
    }
    uint jBVI = 0;
    for(int j = 0; j < pptree.nLevels; j++) {
      const uint bFirst = pptree.tipsVectorIndex[j];
      const uint bLast = pptree.tipsVectorIndex[j + 1] - 1;

      _PRAGMA_OMP_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }
      uint nBranchesDone = 0;
      while(nBranchesDone != bLast - bFirst + 1) {
        const uint unFirst = pptree.branchVectorIndex[jBVI];
        const uint unLast = pptree.branchVectorIndex[jBVI + 1] - 1;

        _PRAGMA_OMP_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.parentNode[i]);
        }

        nBranchesDone +=  pptree.branchVectorIndex[jBVI + 1] - pptree.branchVectorIndex[jBVI];
        ++jBVI;
      }
    }
  }

  void do_pruning_parallel() const {
    spec.initSpecialData();

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < pptree.M - 1; i++) {
    // calculate and store cache information for branch i
    // This is the ideal place to do transformation of branch lengths, etc.
    spec.prepareBranch(i);
  }

  uint jBVI = 0;

  for(int j = 0; j < pptree.nLevels; j++) {
    const uint bFirst = pptree.tipsVectorIndex[j];
    const uint bLast = pptree.tipsVectorIndex[j + 1] - 1;

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }

      uint nBranchesDone = 0;

    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = pptree.branchVectorIndex[jBVI];
      const uint unLast = pptree.branchVectorIndex[jBVI + 1] - 1;

      _PRAGMA_OMP_FOR_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.parentNode[i]);
        }

        nBranchesDone +=  pptree.branchVectorIndex[jBVI + 1] - pptree.branchVectorIndex[jBVI];
      ++jBVI;
#pragma omp barrier
    }
  }
}
  }

  void do_pruning_hybrid() const {
    spec.initSpecialData();

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
    for(uint i = 0; i < pptree.M - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      spec.prepareBranch(i);
    }
  } else {
    // only one (master) thread executes this
    for(uint i = 0; i < pptree.M - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      spec.prepareBranch(i);
    }
  }

  uint jBVI = 0;

  for(int j = 0; j < pptree.nLevels; j++) {
    const uint bFirst = pptree.tipsVectorIndex[j];
    const uint bLast = pptree.tipsVectorIndex[j + 1] - 1;

    if(bLast - bFirst + 1 > nThreads * MIN_CHUNK_SIZE) {
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }
    } else if(tid == 0) {
      // only one (master) thread executes this
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }
    }

    uint nBranchesDone = 0;

    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = pptree.branchVectorIndex[jBVI];
      const uint unLast = pptree.branchVectorIndex[jBVI + 1] - 1;

      if(unLast - unFirst + 1 > nThreads * MIN_CHUNK_SIZE) {
        _PRAGMA_OMP_FOR_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.parentNode[i]);
        }
      } else if(tid == 0) {
        // only one (master) thread executes this
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.parentNode[i]);
        }
      }

      nBranchesDone +=  pptree.branchVectorIndex[jBVI + 1] - pptree.branchVectorIndex[jBVI];
      ++jBVI;
#pragma omp barrier
    }
  }
}
  }
};
}
#endif // ParallelPruning_ParallelPruningAlgorithm_H_
