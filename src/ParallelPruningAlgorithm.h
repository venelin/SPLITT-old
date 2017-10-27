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
//        constructor ParallelPruningMeta(N, branch_starts_, branch_ends_, t).
//       4.3.2 void prepareBranch(uint i)
//       4.3.3 void pruneBranch(uint i)
//       4.3.4 void addToParent(uint i, uint iParent)
//    4.4 If you allocated dynamic resources or memory, you might need to implement
//      a destructor ~MyPruningAlgorithm.
// 5. In your user code, create instances of your class, set the data fields and
//  and call do_pruning() on various parameter sets. After calling do_pruning(),
//  the pruning result at index num_all_nodes_-1 should correspond to the final result at the
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
#include <stack>

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

typedef unsigned int uint;
typedef std::vector<uint> uvec;
typedef std::vector<double> vec;
typedef std::vector<bool> bvec;

// define an NA constant;
const uint NA_UINT = std::numeric_limits<uint>::max();

template <class VectorClass>
inline std::vector<uint> SortIndices(VectorClass const& v) {

  // initialize original index locations
  std::vector<uint> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template<class VectorValues, class VectorPositions>
inline VectorValues At(VectorValues const& v, VectorPositions const& positions) {
  VectorValues sub;
  sub.resize(positions.size());

  size_t sub_i = 0;
  for(auto pit = positions.begin(); pit != positions.end(); pit++,sub_i++){
    sub[sub_i] = v[*pit];
  }
  return sub;
}

template<class VectorValues>
inline VectorValues At(VectorValues const& v, bvec const& mask) {
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
inline uvec Match(uvec const& x, uvec const& table) {
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

inline uvec Seq(uint first, uint last) {
  uvec res(last-first+1);
  std::iota(res.begin(), res.end(), first);
  return res;
}

inline bvec IsNA(uvec const& x) {
  bvec res(x.size(), false);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA_UINT) res[i] = true;
  }
  return res;
}

inline bvec NotIsNA(uvec const& x) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA_UINT) res[i] = false;
  }
  return res;
}

// A tree topology is a tree without branch-lengths. Still, for our purposes,
// it has named nodes. Permuting the nodes results in a different topology.
template<class Node, class BranchLength>
class Tree {
protected:
  uint num_tips_;
  uint num_all_nodes_;
  uvec id_parent_node_;

  uvec id_branch_original_;

  std::map<Node, uint> map_node_to_id_;
  std::vector<Node> map_id_to_node_;
  std::vector<BranchLength> branch_lengths_;

public:
  Tree(std::vector<Node> const& branch_start_nodes,
       std::vector<Node> const& branch_end_nodes,
       std::vector<BranchLength> const& branch_lengths) {

    if(branch_start_nodes.size() != branch_end_nodes.size()) {
      std::ostringstream oss;
      oss<<"branch_start_nodes and branch_end_nodes should be the same size, but were "
         <<branch_start_nodes.size()<<" and "<<
      branch_end_nodes.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    // There should be exactly num_all_nodes_ = number-of-branches + 1
    // distinct nodes
    // This is because each branch can be mapped to its ending node. The +1
    // corresponds to the root node, to which no branch points.
    this->num_all_nodes_ = branch_start_nodes.size() + 1;

    // we distinguish three types of nodes:
    enum NodeType { ROOT, INTERNAL, TIP };

    // initially, we traverse the list of branches and order the nodes as they
    // appear during the traversal from 0 to num_all_nodes_-1. This order is done by
    // incrementing node_id_temp.
    uint node_id_temp = 0;

    std::vector<NodeType> node_types(num_all_nodes_, ROOT);
    this->map_id_to_node_ = std::vector<Node>(num_all_nodes_);
    uvec branch_starts_temp(branch_start_nodes.size(), NA_UINT);
    uvec branch_ends_temp(branch_start_nodes.size(), NA_UINT);
    uvec ending_at(num_all_nodes_ - 1, NA_UINT);

    for(uint i = 0; i < branch_start_nodes.size(); ++i) {
      if(branch_start_nodes[i] == branch_end_nodes[i]) {
        std::ostringstream oss;
        oss<<"Found a branch with the same start and end node ("<<
          branch_start_nodes[i]<<"). Not allowed. ";
        throw std::logic_error(oss.str());
      }

      auto it1 = map_node_to_id_.find(branch_start_nodes[i]);
      if(it1 == map_node_to_id_.end()) {
        // node encountered for the first time
        map_node_to_id_[branch_start_nodes[i]] = node_id_temp; // insert branch_start_nodes[i] in the map
        map_id_to_node_[node_id_temp] = branch_start_nodes[i];
        if(node_types[node_id_temp] == TIP) {
          node_types[node_id_temp] = INTERNAL;
        }
        branch_starts_temp[i] = node_id_temp;
        node_id_temp++;
      } else {
        // node encountered in a previous branch
        if(node_types[it1->second] == TIP) {
          // the previous encounter of the node was as a branch-end
          node_types[it1->second] = INTERNAL;
        } else {
          // do nothing
        }
        branch_starts_temp[i] = it1->second;
      }

      auto it2 = map_node_to_id_.find(branch_end_nodes[i]);
      if(it2 == map_node_to_id_.end()) {
        // node encountered for the first time
        map_node_to_id_[branch_end_nodes[i]] = node_id_temp;
        map_id_to_node_[node_id_temp] = branch_end_nodes[i];

        if(node_types[node_id_temp] == ROOT) {
          // not known if the node has descendants, so we set its type to TIP.
          node_types[node_id_temp] = TIP;
        }
        branch_ends_temp[i] = node_id_temp;
        ending_at[node_id_temp] = i;
        node_id_temp++;
      } else {
        // node has been previously encountered
        if(ending_at[it2->second] != NA_UINT) {
          std::ostringstream oss;
          oss<<"Found at least two branches ending at the same node ("<<
            it2->first<<"). Check for cycles or repeated branches. ";
          throw std::logic_error(oss.str());
        } else {
          if(node_types[it2->second] == ROOT) {
            // the previous enounters of the node were as branch-start -> set
            // the node's type to INTERNAL, because we know for sure that it
            // has descendants.
            node_types[it2->second] = INTERNAL;
          }
          branch_ends_temp[i] = it2->second;
          ending_at[it2->second] = i;
        }
      }
    }

    if(map_node_to_id_.size() != num_all_nodes_) {
      std::ostringstream oss;
      oss<<"The number of distinct nodes ("<<map_node_to_id_.size()<<
        ") should equal the number-of-branches+1 ("<<num_all_nodes_<<").";
      throw std::logic_error(oss.str());
    }

    auto num_roots = count(node_types.begin(), node_types.end(), ROOT);
    if(num_roots != 1) {
      std::ostringstream oss;
      oss<<"There should be exactly one ROOT node, but "<<num_roots<<
        " were found. Check for cycles or for multiple trees.";
      throw std::logic_error(oss.str());
    }

    this->num_tips_ = count(node_types.begin(), node_types.end(), TIP);
    if(num_tips_ == 0) {
      std::ostringstream oss;
      oss<<"There should be at least one TIP node, but none"<<
        " was found. Check for cycles.";
      throw std::logic_error(oss.str());
    }

    // assign new ids according to the following convention:
    // tips are numbered from 0 to num_tips_ - 1;
    // internal nodes are numbered from num_tips_ to num_all_nodes_ - 2;
    // root is numbered num_all_nodes_ - 1;
    std::vector<uint> node_ids(num_all_nodes_, NA_UINT);

    uint tip_no = 0, internal_no = num_tips_;
    for(uint i = 0; i < num_all_nodes_; i++) {
      if(node_types[i] == TIP) {
        node_ids[i] = tip_no;
        tip_no ++;
      } else if(node_types[i] == INTERNAL) {
        node_ids[i] = internal_no;
        internal_no++;
      } else {
        // Here node_types[i] == ROOT should be true
        node_ids[i] = num_all_nodes_ - 1;
      }
      map_node_to_id_[map_id_to_node_[i]] = node_ids[i];
    }

    this->map_id_to_node_ = At(map_id_to_node_, SortIndices(node_ids));

    this->id_parent_node_ = uvec(num_all_nodes_ - 1);
    this->branch_lengths_ = std::vector<BranchLength>(num_all_nodes_ - 1);

    this->id_branch_original_ = uvec(num_all_nodes_ - 1);

    for(uint i = 0; i < num_all_nodes_ - 1; i++) {
      uint branch_start_i = node_ids[branch_starts_temp[i]];
      uint branch_end_i = node_ids[branch_ends_temp[i]];
      id_parent_node_[branch_end_i] = branch_start_i;
      branch_lengths_[branch_end_i] = branch_lengths[i];
      id_branch_original_[branch_end_i] = i;
    }
  }

  uint num_all_nodes() const {
    return num_all_nodes_;
  }

  uint num_tips() const {
    return num_tips_;
  }

  BranchLength branch_length(uint i) const {
    return branch_lengths_[i];
  }

  // return the node at position id;
  Node const& node(uint id) const {
    return map_id_to_node_[id];
  }

  // get the internally stored id of node
  uint id_node(Node const& node) const {
    auto it = map_node_to_id_.find(node);
    if(it == map_node_to_id_.end()) {
      return NA_UINT;
    } else {
      return it->second;
    }
  }
};

template<class Node, class BranchLength>
class ParallelPruningTree: public Tree<Node, BranchLength> {
protected:
  // number of pruning steps
  uint num_levels_;

  uvec nodes_to_prune_;
  uvec nodes_to_update_parent_;

  uvec orderNodeIds;

  void CreatePruningOrder() {
    // insert a fictive branch leading to the root of the tree.
    uvec branch_ends = Seq(0, this->num_all_nodes_ - 1);

    uvec num_children_remaining(this->num_all_nodes_, 0);
    for(uint i : this->id_parent_node_) num_children_remaining[i]++;

    this->nodes_to_prune_ = uvec(1, 0);
    this->nodes_to_update_parent_ = uvec(1, 0);

    // start by pruning the tips of the tree
    uvec tips_this_level = Seq(0, this->num_tips_ - 1);

    uvec default_pos_vector;
    default_pos_vector.reserve(2);
    std::vector<uvec> pos_of_parent(this->num_all_nodes_ - this->num_tips_,
                                    default_pos_vector);

    while(tips_this_level[0] != this->num_all_nodes_ - 1) { // while the root has not become a tip itself

      nodes_to_prune_.push_back(
        nodes_to_prune_[nodes_to_prune_.size() - 1] + tips_this_level.size());


      // unique parents at this level
      uvec parents_this_level;
      parents_this_level.reserve(tips_this_level.size());

      // empty the tips_this_level vector so it can be filled in with new tips to be pruned
      //tips_this_level.clear();
      uvec tips_next_level;

      for(uint i = 0; i < tips_this_level.size(); i++) {
        uint i_parent = this->id_parent_node_[tips_this_level[i]];
        if(pos_of_parent[i_parent - this->num_tips_].empty()) {
          parents_this_level.push_back(i_parent);
        }
        pos_of_parent[i_parent - this->num_tips_].push_back(i);
      }

      uint num_parents_remaining = parents_this_level.size();
      while( num_parents_remaining ) {

        uint num_parent_updates = 0;
        for(auto i_parent: parents_this_level) {
          if(!pos_of_parent[i_parent - this->num_tips_].empty()) {
            uint i = pos_of_parent[i_parent - this->num_tips_].back();

            num_parent_updates ++;
            this->orderNodeIds.push_back(tips_this_level[i]);

            num_children_remaining[i_parent]--;
            if(num_children_remaining[i_parent] == 0) {
              tips_next_level.push_back(i_parent);
            }
            pos_of_parent[i_parent - this->num_tips_].pop_back();
            if(pos_of_parent[i_parent - this->num_tips_].empty()) {
              num_parents_remaining--;
            }
          }
        }

        nodes_to_update_parent_.push_back(
          nodes_to_update_parent_[nodes_to_update_parent_.size() - 1] + num_parent_updates);
      }

      tips_this_level = tips_next_level;
    }

    this->num_levels_ = nodes_to_prune_.size() - 1;
    this->orderNodeIds.push_back(this->num_all_nodes_ - 1);

    auto orderBranches = get_orderBranches();
    this->branch_lengths_ = At(this->branch_lengths_, orderBranches);

    uvec branchEndsReord = orderBranches;
    branchEndsReord.push_back(this->num_all_nodes_ - 1);
    this->id_parent_node_ = Match(At(this->id_parent_node_, orderBranches),
                               branchEndsReord);
  }

  // private default constructor;
  ParallelPruningTree() {};


public:
  ParallelPruningTree(
    std::vector<Node> const& branch_start_nodes,
    std::vector<Node> const& branch_end_nodes,
    std::vector<BranchLength> const& branch_lengths): Tree<Node, BranchLength>(branch_start_nodes, branch_end_nodes, branch_lengths) {

    CreatePruningOrder();
  }

  uint num_levels() const {
    return num_levels_;
  }

  uint ParentPruneIndex(uint child_prune_index) const {
    return this->id_parent_node_[child_prune_index];
  }

  uvec const& nodes_to_prune() const {
    return nodes_to_prune_;
  }

  uvec const& nodes_to_update_parent() const {
    return nodes_to_update_parent_;
  }

  // The order in which the node ids get processed during the pruning, including
  // tips, internal and root nodes. Each element of the returned
  // vector is a node-id. See also get_orderNodes.
  uvec get_orderNodeIds() const {
    return orderNodeIds;
  }

  // The order in which the branches get processed during the pruning.
  uvec get_orderBranches() const {
    uvec orderBranchEnds(orderNodeIds.begin(), orderNodeIds.end() - 1);
    return Match(orderBranchEnds, Seq(0, this->num_all_nodes_ - 1));
  }

  // A root-to-node distance vector in the order of pruning processing
  vec get_nodeHeights() const {
    vec h(this->num_all_nodes_, 0);
    for(int i = this->num_all_nodes_ - 2; i >= 0; i--) {
      h[i] = h[this->id_parent_node_[i]] + this->branch_lengths_[i];
    }
    return h;
  }

  // returns a vector of positions in nodes in the order of pruning.
  // See also get_orderNodeIds.
  uvec order_nodes(std::vector<Node> const& nodes) const {
    uvec ids(nodes.size());
    for(uint i = 0; i < nodes.size(); ++i) {
      auto it = this->map_node_to_id_.find(nodes[i]);
      if(it == this->map_node_to_id_.end()) {
        std::ostringstream oss;
        oss<<"At least one of the nodes is not present in the tree ("<<
          nodes[i]<<").";
        throw std::invalid_argument(oss.str());
      } else {
        ids[i] = it->second;
      }
    }
    uvec m = Match(orderNodeIds, ids);
    return At(m, NotIsNA(m));
  }
};

//  a skeleton for a ParallelPruningSpec class
// Curiously recurring template pattern -  static polymorphism
template<class ParallelPruningTree, class PruningSpec>
class ParallelPruningAlgorithm {
protected:
  uint nThreads;
  const ParallelPruningTree& pptree;
  PruningSpec& spec;

  // a for loop in do_pruning_hybrid will be omp-for-ed if the number of
  // iterations is bigger than nThreads * minChunkSizeForHybrid.
  uint minChunkSizeForHybrid = 2;
  uint bestMinChunkSizeForHybrid = 2;
  bool minChunkSizeGrowing = 1;
  enum PruningModeAuto { SINGLE_THREAD = 1, MULTI_THREAD = 2, HYBRID = 3 };

  uint pmaCurrentStep = 0;
  const uint PRUNING_MODE_AUTO_MAX_STEPS = 20;

  double pmaMinDuration = std::numeric_limits<double>::max();
  double pmaHybridMinDuration = std::numeric_limits<double>::max();
  std::vector<double> pmaTuningDurations;
  PruningModeAuto pmaBestMode = MULTI_THREAD;
  std::vector<PruningModeAuto> pmaTuningModes;
  std::vector<uint> pmaTuningMinChunkSizes;
public:
  ParallelPruningAlgorithm(const ParallelPruningTree& _pptree, PruningSpec& _spec):
  pptree(_pptree), spec(_spec),
  pmaTuningModes(PRUNING_MODE_AUTO_MAX_STEPS + 1, HYBRID),
  pmaTuningDurations(
    PRUNING_MODE_AUTO_MAX_STEPS + 1, std::numeric_limits<double>::max()),
  pmaTuningMinChunkSizes(PRUNING_MODE_AUTO_MAX_STEPS + 1, 2) {

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
    minChunkSizeForHybrid = 2;
    pmaCurrentStep = 0;
    pmaTuningModes[0] = SINGLE_THREAD;
    pmaTuningModes[1] = MULTI_THREAD;
    pmaTuningModes[11] = SINGLE_THREAD;
    pmaTuningModes[PRUNING_MODE_AUTO_MAX_STEPS-1] = SINGLE_THREAD;
  }

  uint get_nThreads() const {
    return nThreads;
  }

  uint get_bestMinChunkSizeForHybrid() const {
    return bestMinChunkSizeForHybrid;
  }


  std::vector<double>  get_pmaTuningDurations() const {
    return pmaTuningDurations;
  }

  PruningModeAuto get_pmaBestMode() const {
    return pmaBestMode;
  }

  std::vector<uint> get_pmaTuningMinChunkSizes() const {
    return pmaTuningMinChunkSizes;
  }

  void do_pruning(int mode) {
    switch(mode) {
    case 0: do_pruning_auto(); break;
    case 1: do_pruning_single_thread(); break;
    case 2: do_pruning_multi_thread(); break;
    case 3: do_pruning_hybrid(); break;
    default: do_pruning_auto();
    }
  }

protected:
  void do_pruning_auto() {
    PruningModeAuto mode;
    std::chrono::steady_clock::time_point start, end;
    double duration;
    bool isTuning = false;
    if(pmaCurrentStep < PRUNING_MODE_AUTO_MAX_STEPS) {
      mode = pmaTuningModes[pmaCurrentStep];
      isTuning = true;
    } else {
      mode = pmaBestMode;
      minChunkSizeForHybrid = bestMinChunkSizeForHybrid;
      isTuning = false;
    }

    if( isTuning ) {
      switch(mode) {

      case SINGLE_THREAD:
        start = std::chrono::steady_clock::now();
        do_pruning(1);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double, std::milli>(end - start).count();
        pmaTuningDurations[pmaCurrentStep] = duration;
        if(duration < pmaMinDuration) {
          pmaMinDuration = duration;
          pmaBestMode = SINGLE_THREAD;
        }
        pmaCurrentStep++;
        break;

      case MULTI_THREAD:
        start = std::chrono::steady_clock::now();
        do_pruning(2);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double, std::milli>(end - start).count();
        pmaTuningDurations[pmaCurrentStep] = duration;
        if(duration < pmaMinDuration) {
          pmaMinDuration = duration;
          pmaBestMode = MULTI_THREAD;
        }
        pmaCurrentStep++;
        break;

      case HYBRID:
        start = std::chrono::steady_clock::now();
        do_pruning(3);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double, std::milli>(end - start).count();
        pmaTuningDurations[pmaCurrentStep] = duration;
        pmaTuningMinChunkSizes[pmaCurrentStep] = minChunkSizeForHybrid;

        if(duration < pmaMinDuration) {
          pmaMinDuration = duration;
          pmaBestMode = HYBRID;
        }
        if(duration < pmaHybridMinDuration) {
          pmaHybridMinDuration = duration;
          bestMinChunkSizeForHybrid = minChunkSizeForHybrid;
        } else if(duration > pmaTuningDurations[pmaCurrentStep - 1]) {
          // not improving from the current step, so change direction of minChunkSize
          if(!minChunkSizeGrowing) {
            minChunkSizeGrowing = 1;
          } else if(minChunkSizeGrowing & minChunkSizeForHybrid > 2) {
            // if minChunkSizeForHybrid is hitting the bottom, keep growing
            minChunkSizeGrowing = 0;
          }
        }
        minChunkSizeForHybrid += minChunkSizeForHybrid * minChunkSizeGrowing -
          (minChunkSizeForHybrid / 4) * (!minChunkSizeGrowing);
        pmaCurrentStep++;
        break;
      }
    } else {
      switch(mode) {
      case HYBRID: do_pruning(3); pmaCurrentStep++; break;
      case SINGLE_THREAD: do_pruning(1); pmaCurrentStep++; break;
      case MULTI_THREAD: do_pruning(2); pmaCurrentStep++; break;
      }
    }

  }

  void do_pruning_single_thread() {
    spec.initSpecialData();

    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < pptree.num_all_nodes() - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      spec.prepareBranch(i);
    }

    uint jBI = 0;
    for(int j = 0; j < pptree.num_levels(); j++) {
      const uint bFirst = pptree.nodes_to_prune()[j];
      const uint bLast = pptree.nodes_to_prune()[j + 1] - 1;

      _PRAGMA_OMP_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }
      uint nBranchesDone = 0;
      while(nBranchesDone != bLast - bFirst + 1) {
        const uint unFirst = pptree.nodes_to_update_parent()[jBI];
        const uint unLast = pptree.nodes_to_update_parent()[jBI + 1] - 1;

        _PRAGMA_OMP_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.ParentPruneIndex(i));
        }

        nBranchesDone +=  pptree.nodes_to_update_parent()[jBI + 1] - pptree.nodes_to_update_parent()[jBI];
        ++jBI;
      }
    }
  }

  void do_pruning_multi_thread() {
    spec.initSpecialData();

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < pptree.num_all_nodes() - 1; i++) {
    // calculate and store cache information for branch i
    // This is the ideal place to do transformation of branch lengths, etc.
    spec.prepareBranch(i);
  }

  uint jBI = 0;

  for(int j = 0; j < pptree.num_levels(); j++) {
    const uint bFirst = pptree.nodes_to_prune()[j];
    const uint bLast = pptree.nodes_to_prune()[j + 1] - 1;

    _PRAGMA_OMP_FOR_SIMD
    for(uint i = bFirst; i < bLast + 1; i++) {
      // perform the main calculation for branch i based on the pruning
      // results from its daughter branches.
      spec.pruneBranch(i);
    }

    uint nBranchesDone = 0;

    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = pptree.nodes_to_update_parent()[jBI];
      const uint unLast = pptree.nodes_to_update_parent()[jBI + 1] - 1;

      _PRAGMA_OMP_FOR_SIMD
      for(uint i = unFirst; i < unLast + 1; i++) {
        // store or add up the result from branch i to the results from its
        // sibling branches, so that these results can be used for the
        // pruning of the parent branch.
        spec.addToParent(i, pptree.ParentPruneIndex(i));
      }

      nBranchesDone +=  pptree.nodes_to_update_parent()[jBI + 1] - pptree.nodes_to_update_parent()[jBI];
      ++jBI;
#pragma omp barrier
    }
  }
}
  }

  void do_pruning_hybrid() {
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
    for(uint i = 0; i < pptree.num_all_nodes() - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      spec.prepareBranch(i);
    }
  } else {
    // only one (master) thread executes this
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < pptree.num_all_nodes() - 1; i++) {
      // calculate and store cache information for branch i
      // This is the ideal place to do transformation of branch lengths, etc.
      spec.prepareBranch(i);
    }
  }

  uint jBI = 0;

  for(int j = 0; j < pptree.num_levels(); j++) {
    const uint bFirst = pptree.nodes_to_prune()[j];
    const uint bLast = pptree.nodes_to_prune()[j + 1] - 1;

    if(bLast - bFirst + 1 > nThreads * minChunkSizeForHybrid) {
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }
    } else if(tid == 0) {
      // only one (master) thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = bFirst; i < bLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // results from its daughter branches.
        spec.pruneBranch(i);
      }
    }

    uint nBranchesDone = 0;

    while(nBranchesDone != bLast - bFirst + 1) {
      const uint unFirst = pptree.nodes_to_update_parent()[jBI];
      const uint unLast = pptree.nodes_to_update_parent()[jBI + 1] - 1;

      if(unLast - unFirst + 1 > nThreads * minChunkSizeForHybrid) {
        _PRAGMA_OMP_FOR_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.ParentPruneIndex(i));
        }
      } else if(tid == 0) {
        // only one (master) thread executes this
        _PRAGMA_OMP_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          // store or add up the result from branch i to the results from its
          // sibling branches, so that these results can be used for the
          // pruning of the parent branch.
          spec.addToParent(i, pptree.ParentPruneIndex(i));
        }
      }

      nBranchesDone +=  pptree.nodes_to_update_parent()[jBI + 1] - pptree.nodes_to_update_parent()[jBI];
      ++jBI;
#pragma omp barrier
    }
  }
}
  }
};
}
#endif // ParallelPruning_ParallelPruningAlgorithm_H_
