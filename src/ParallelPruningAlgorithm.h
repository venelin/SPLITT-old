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
//  and call DoPruning() on various parameter sets. After calling DoPruning(),
//  the pruning result at index num_nodes_-1 should correspond to the final result at the
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
#include <mutex>

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
template<class Node, class Length>
class Tree {
protected:
  uint num_tips_;
  uint num_nodes_;
  uvec id_parent_;

  uvec id_branch_original_;

  std::map<Node, uint> map_node_to_id_;
  std::vector<Node> map_id_to_node_;
  std::vector<Length> lengths_;

public:
  Tree(std::vector<Node> const& branch_start_nodes,
       std::vector<Node> const& branch_end_nodes,
       std::vector<Length> const& branch_lengths) {

    if(branch_start_nodes.size() != branch_end_nodes.size()) {
      std::ostringstream oss;
      oss<<"branch_start_nodes and branch_end_nodes should be the same size, but were "
         <<branch_start_nodes.size()<<" and "<<
      branch_end_nodes.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    // There should be exactly num_nodes_ = number-of-branches + 1
    // distinct nodes
    // This is because each branch can be mapped to its ending node. The +1
    // corresponds to the root node, to which no branch points.
    this->num_nodes_ = branch_start_nodes.size() + 1;

    // we distinguish three types of nodes:
    enum NodeType { ROOT, INTERNAL, TIP };

    // initially, we traverse the list of branches and order the nodes as they
    // appear during the traversal from 0 to num_nodes_-1. This order is done by
    // incrementing node_id_temp.
    uint node_id_temp = 0;

    std::vector<NodeType> node_types(num_nodes_, ROOT);
    this->map_id_to_node_ = std::vector<Node>(num_nodes_);
    uvec branch_starts_temp(branch_start_nodes.size(), NA_UINT);
    uvec branch_ends_temp(branch_start_nodes.size(), NA_UINT);
    uvec ending_at(num_nodes_ - 1, NA_UINT);

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

    if(map_node_to_id_.size() != num_nodes_) {
      std::ostringstream oss;
      oss<<"The number of distinct nodes ("<<map_node_to_id_.size()<<
        ") should equal the number-of-branches+1 ("<<num_nodes_<<").";
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
    // internal nodes are numbered from num_tips_ to num_nodes_ - 2;
    // root is numbered num_nodes_ - 1;
    std::vector<uint> node_ids(num_nodes_, NA_UINT);

    uint tip_no = 0, internal_no = num_tips_;
    for(uint i = 0; i < num_nodes_; i++) {
      if(node_types[i] == TIP) {
        node_ids[i] = tip_no;
        tip_no ++;
      } else if(node_types[i] == INTERNAL) {
        node_ids[i] = internal_no;
        internal_no++;
      } else {
        // Here node_types[i] == ROOT should be true
        node_ids[i] = num_nodes_ - 1;
      }
      map_node_to_id_[map_id_to_node_[i]] = node_ids[i];
    }

    this->map_id_to_node_ = At(map_id_to_node_, SortIndices(node_ids));

    this->id_parent_ = uvec(num_nodes_ - 1);
    this->lengths_ = std::vector<Length>(num_nodes_ - 1);

    this->id_branch_original_ = uvec(num_nodes_ - 1);

    for(uint i = 0; i < num_nodes_ - 1; i++) {
      uint branch_start_i = node_ids[branch_starts_temp[i]];
      uint branch_end_i = node_ids[branch_ends_temp[i]];
      id_parent_[branch_end_i] = branch_start_i;
      lengths_[branch_end_i] = branch_lengths[i];
      id_branch_original_[branch_end_i] = i;
    }
  }

  uint num_nodes() const {
    return num_nodes_;
  }

  uint num_tips() const {
    return num_tips_;
  }

  Length LengthOfBranch(uint i) const {
    return lengths_[i];
  }

  // return the node at position id;
  Node const& FindNodeWithId(uint id) const {
    return map_id_to_node_[id];
  }

  // get the internally stored id of node
  uint FindIdOfNode(Node const& node) const {
    auto it = map_node_to_id_.find(node);
    if(it == map_node_to_id_.end()) {
      return NA_UINT;
    } else {
      return it->second;
    }
  }

  uint FindIdOfParent(uint id_child) const {
    return this->id_parent_[id_child];
  }
};

template<class Node, class Length>
class ParallelPruningTree: public Tree<Node, Length> {
protected:
  uvec ranges_id_prune_;
  uvec ranges_id_update_parent_;

  // private default constructor;
  ParallelPruningTree() {}
public:
  ParallelPruningTree(
    std::vector<Node> const& branch_start_nodes,
    std::vector<Node> const& branch_end_nodes,
    std::vector<Length> const& branch_lengths): Tree<Node, Length>(branch_start_nodes, branch_end_nodes, branch_lengths) {

      // insert a fictive branch leading to the root of the tree.
      uvec branch_ends = Seq(0, this->num_nodes_ - 1);

      uvec num_children_remaining(this->num_nodes_, 0);
      for(uint i : this->id_parent_) num_children_remaining[i]++;

      this->ranges_id_prune_ = uvec(1, 0);
      this->ranges_id_update_parent_ = uvec(1, 0);

      // start by pruning the tips of the tree
      uvec tips_this_level = Seq(0, this->num_tips_ - 1);

      uvec default_pos_vector;
      default_pos_vector.reserve(2);
      std::vector<uvec> pos_of_parent(this->num_nodes_ - this->num_tips_,
                                      default_pos_vector);

      uvec order_branches;
      order_branches.reserve(this->num_nodes_);

      while(tips_this_level[0] != this->num_nodes_ - 1) {
        // while the root has not become a tip itself
        ranges_id_prune_.push_back(
          ranges_id_prune_[ranges_id_prune_.size() - 1] + tips_this_level.size());


        // unique parents at this level
        uvec parents_this_level;
        parents_this_level.reserve(tips_this_level.size());

        uvec tips_next_level;
        tips_next_level.reserve(tips_this_level.size() / 2);

        for(uint i = 0; i < tips_this_level.size(); i++) {
          uint i_parent = this->id_parent_[tips_this_level[i]];
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
              order_branches.push_back(tips_this_level[i]);

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

          ranges_id_update_parent_.push_back(
            ranges_id_update_parent_[ranges_id_update_parent_.size() - 1] + num_parent_updates);
        }

        tips_this_level = tips_next_level;
      }

      this->lengths_ = At(this->lengths_, order_branches);

      uvec id_old = order_branches;
      id_old.push_back(this->num_nodes_ - 1);

      this->id_parent_ = Match(At(this->id_parent_, order_branches),
                               id_old);

      // update maps
      std::vector<Node> map_id_to_node(this->num_nodes_);
      for (uint i = 0; i < this->num_nodes_; i++) {
        map_id_to_node[i] = this->map_id_to_node_[id_old[i]];
        this->map_node_to_id_[map_id_to_node[i]] = i;
      }
      std::swap(this->map_id_to_node_, map_id_to_node);
    }

  uint num_levels() const {
    return ranges_id_prune_.size() - 1;
  }

  uvec const& ranges_id_prune() const {
    return ranges_id_prune_;
  }

  std::pair<uint, uint> RangeIdPrune(uint i_level) const {
    return std::pair<uint, uint>(ranges_id_prune_[i_level],
                                 ranges_id_prune_[i_level+1] - 1);
  }

  uvec const& ranges_id_update_parent() const {
    return ranges_id_update_parent_;
  }

  std::pair<uint, uint> RangeIdUpdateParent(uint i_step) const {
    return std::pair<uint, uint>(ranges_id_update_parent_[i_step],
                                 ranges_id_update_parent_[i_step+1] - 1);
  }

  // A root-to-node distance vector in the order of pruning processing
  std::vector<Length> CalculateHeights() const {
    std::vector<Length> h(this->num_nodes_, 0);
    for(int i = this->num_nodes_ - 2; i >= 0; i--) {
      h[i] = h[this->id_parent_[i]] + this->lengths_[i];
    }
    return h;
  }

  // returns a vector of positions in nodes in the order of pruning.
  // See also get_orderNodeIds.
  uvec OrderNodes(std::vector<Node> const& nodes) const {
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
    uvec m = Match(Seq(0, this->num_nodes_ - 1), ids);
    return At(m, NotIsNA(m));
  }
};

//  a skeleton for a ParallelPruningSpec class
// Curiously recurring template pattern -  static polymorphism
template<class ParallelPruningTree, class PruningSpec>
class ParallelPruningAlgorithm {
protected:
  const ParallelPruningTree& pptree_;
  PruningSpec& spec_;


  uint num_threads_;

  uint current_step_tuning_ = 0;
  uint fastest_step_tuning_ = 0;

  double min_duration_tuning_ = std::numeric_limits<double>::max();
  std::vector<double> durations_tuning_;

  enum PruningModeAuto { SINGLE_THREAD = 1, MULTI_THREAD = 2, HYBRID = 3 };

  const uvec min_sizes_chunk_ = {1, 2, 4, 8, 16, 32, 64};

public:
  ParallelPruningAlgorithm(
    const ParallelPruningTree& pptree, PruningSpec& spec):
  pptree_(pptree), spec_(spec) {

#ifdef _OPENMP
#pragma omp parallel
{
  uint tid = omp_get_thread_num();
  // only master thread does this
  if(tid == 0) {
    this->num_threads_ = omp_get_num_threads();
  }
}
#else
this->num_threads_ = 1;
#endif // #ifdef _OPENMP
  }

  uint num_threads() const {
    return num_threads_;
  }

  bool IsTuning() const {
    return current_step_tuning_ < min_sizes_chunk_.size() * min_sizes_chunk_.size() + 4;
  }

  uint IndexMinSizeChunkPrune() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ((step - 4) / min_sizes_chunk_.size()) % min_sizes_chunk_.size();
  }

  uint IndexMinSizeChunkUpdate() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ((step - 4) % min_sizes_chunk_.size()) % min_sizes_chunk_.size();
  }

  PruningModeAuto ModeTuning() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    if(step < 2) {
      return SINGLE_THREAD;
    } else if(step < 4) {
      return MULTI_THREAD;
    } else {
      return HYBRID;
    }
  }

  uint min_size_chunk_prune() const {
    return min_sizes_chunk_[IndexMinSizeChunkPrune()];
  }

  uint min_size_chunk_update() const {
    return min_sizes_chunk_[IndexMinSizeChunkUpdate()];
  }

  uint fastest_step_tuning() const {
    return fastest_step_tuning_;
  }

  std::vector<double>  durations_tuning() const {
    return durations_tuning_;
  }

  void DoPruning(int mode) {
    switch(mode) {
    case 0: DoPruningAuto(); break;
    case 1: DoPruningSingleThread(); break;
    case 2: DoPruningMultiThread(); break;
    case 3: DoPruningHybrid(); break;
    default: DoPruningAuto();
    }
  }

protected:
  void DoPruningAuto() {
    PruningModeAuto mode;
    std::chrono::steady_clock::time_point start, end;
    double duration;

    if( /*IsTuning()*/ false ) {
      auto mode = ModeTuning();

      switch(mode) {
      case SINGLE_THREAD:
        start = std::chrono::steady_clock::now();
        DoPruning(1);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double, std::milli>(end - start).count();
        break;

      case MULTI_THREAD:
        start = std::chrono::steady_clock::now();
        DoPruning(2);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double, std::milli>(end - start).count();
        break;

      case HYBRID:
        start = std::chrono::steady_clock::now();
        DoPruning(3);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double, std::milli>(end - start).count();
        break;
      }

      durations_tuning_.push_back(duration);
      if(duration < min_duration_tuning_) {
        min_duration_tuning_ = duration;
        fastest_step_tuning_ = current_step_tuning_;
      }
      current_step_tuning_++;

    } else {
      auto mode = HYBRID; //ModeTuning();

      switch(mode) {
      case HYBRID: DoPruning(3); break;
      case SINGLE_THREAD: DoPruning(1); break;
      case MULTI_THREAD: DoPruning(2); break;
      }
    }

  }

  void DoPruningSingleThread() {
    spec_.initSpecialData();

    _PRAGMA_OMP_SIMD
      for(uint i = 0; i < pptree_.num_nodes() - 1; i++) {
        spec_.prepareBranch(i);
      }

      uint i_update = 0;
    for(uint i_level = 0; i_level < pptree_.num_levels(); i_level++) {
      auto range_prune = pptree_.RangeIdPrune(i_level);
      _PRAGMA_OMP_SIMD
        for(uint i = range_prune.first; i <= range_prune.second; i++) {
          spec_.pruneBranch(i);
        }

        uint num_branches_done = 0;
      while(num_branches_done != range_prune.second - range_prune.first + 1) {
        auto range_update = pptree_.RangeIdUpdateParent(i_update);

        _PRAGMA_OMP_SIMD
          for(uint i = range_update.first; i <= range_update.second; i++) {
            spec_.addToParent(i, pptree_.FindIdOfParent(i));
          }

          num_branches_done +=  range_update.second - range_update.first + 1;
        ++i_update;
      }
    }
  }

  void DoPruningMultiThread() {
    spec_.initSpecialData();

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < pptree_.num_nodes() - 1; i++) {
    spec_.prepareBranch(i);
  }

  uint i_update = 0;
  for(uint i_level = 0; i_level < pptree_.num_levels(); i_level++) {
    auto range_prune = pptree_.RangeIdPrune(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_prune.first; i <= range_prune.second; i++) {
        spec_.pruneBranch(i);
      }

      uint num_branches_done = 0;
    while(num_branches_done != range_prune.second - range_prune.first + 1) {
      auto range_update = pptree_.RangeIdUpdateParent(i_update);

      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_update.first; i <= range_update.second; i++) {
          spec_.addToParent(i, pptree_.FindIdOfParent(i));
        }

        num_branches_done +=  range_update.second - range_update.first + 1;
      ++i_update;
    }
  }
}
  }

  void DoPruningHybrid() {
    spec_.initSpecialData();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif


  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < pptree_.num_nodes() - 1; i++) {
      spec_.prepareBranch(i);
    }

    uint i_update = 0;

  for(int i_level = 0; i_level < pptree_.num_levels(); i_level++) {
    auto range_prune = pptree_.RangeIdPrune(i_level);

    if(range_prune.second - range_prune.first + 1 >
         //num_threads_ * min_size_chunk_prune()) {
         num_threads_ * 32) {

      _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_prune.first; i <= range_prune.second; i++) {
        spec_.pruneBranch(i);
      }

    } else if(tid == 0) {

      // only one (master) thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_prune.first; i <= range_prune.second; i++) {
        spec_.pruneBranch(i);
      }

    }

    uint num_branches_done = 0;

    while(num_branches_done != range_prune.second - range_prune.first + 1) {

      auto range_update = pptree_.RangeIdUpdateParent(i_update);

      if (range_update.second - range_update.first + 1 >
            num_threads_ * 0) {
        //num_threads_ * min_size_chunk_prune()) {

        _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_update.first; i <= range_update.second; i++) {
          spec_.addToParent(i, pptree_.FindIdOfParent(i));
        }

      } else if (tid == 0) {

        // only one (master) thread executes this
        _PRAGMA_OMP_SIMD
        for(uint i = range_update.first; i <= range_update.second; i++) {
          spec_.addToParent(i, pptree_.FindIdOfParent(i));
        }

      }

      num_branches_done +=  range_update.second - range_update.first + 1;
      ++i_update;
#pragma omp barrier
    }
  }
}
  }
};
}
#endif // ParallelPruning_ParallelPruningAlgorithm_H_
