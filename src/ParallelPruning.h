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


#ifndef ParallelPruning_ParallelPruning_H_
#define ParallelPruning_ParallelPruning_H_

#include <algorithm>
#include <vector>
#include <math.h>
#include <sstream>
#include <limits>
#include <numeric>
#include <chrono>
#include <unordered_map>
#include <mutex>
#include <iostream>

#ifdef _OPENMP

// Need to decide wheter to use '#pragma omp for' or '#pragma omp for simd'
#if _OPENMP >= 201307  // OMP 4.0 or higher

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#else // #if _OPENMP >= 201307

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD  /*_Pragma("omp simd")*/

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

const uvec EMPTY_UVEC;

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

template<class Node, class Length>
class Tree {
public:
  typedef Node NodeType;
  typedef Length LengthType;

private:
  // private default constructor
  Tree() {}

protected:
  uint num_tips_;
  uint num_nodes_;
  uvec id_parent_;

  typedef std::unordered_map<NodeType, uint> MapType;
  MapType map_node_to_id_;
  std::vector<NodeType> map_id_to_node_;
  std::vector<LengthType> lengths_;
  std::vector<uvec> id_child_nodes_;

  void init_id_child_nodes() {
    id_child_nodes_ = std::vector<uvec>(this->num_nodes() - this->num_tips());

    // fill child vectors
    for(uint i = 0; i < this->num_nodes() - 1; i++) {
      id_child_nodes_[this->FindIdOfParent(i) - this->num_tips()].push_back(i);
    }
  }

public:
  // pass an empty vector for branch_lengths for a tree without branch lengths.
  Tree(std::vector<NodeType> const& branch_start_nodes,
       std::vector<NodeType> const& branch_end_nodes,
       std::vector<LengthType> const& branch_lengths) {

    std::cout<<"Tree1"<<std::endl;

    if(branch_start_nodes.size() != branch_end_nodes.size()) {
      std::ostringstream oss;
      oss<<"branch_start_nodes and branch_end_nodes should be the same size, but were "
         <<branch_start_nodes.size()<<" and "<<
      branch_end_nodes.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    std::cout<<"Tree2"<<std::endl;
    // There should be exactly num_nodes_ = number-of-branches + 1
    // distinct nodes
    // This is because each branch can be mapped to its ending node. The +1
    // corresponds to the root node, to which no branch points.
    this->num_nodes_ = branch_start_nodes.size() + 1;

    std::cout<<"Tree3"<<std::endl;

    // we distinguish three types of nodes:
    enum NodeRole { ROOT, INTERNAL, TIP };

    // initially, we traverse the list of branches and order the nodes as they
    // appear during the traversal from 0 to num_nodes_-1. This order is done by
    // incrementing node_id_temp.
    uint node_id_temp = 0;

    std::vector<NodeRole> node_types(num_nodes_, ROOT);

    this->map_id_to_node_.resize(num_nodes_);

    this->map_node_to_id_.reserve(num_nodes_);

    uvec branch_starts_temp(branch_start_nodes.size(), NA_UINT);
    uvec branch_ends_temp(branch_start_nodes.size(), NA_UINT);
    uvec ending_at(num_nodes_ - 1, NA_UINT);

    std::vector<typename MapType::iterator> it_map_node_to_id_;
    it_map_node_to_id_.reserve(num_nodes_);

    std::cout<<"Tree4"<<std::endl;

    for(uint i = 0; i < branch_start_nodes.size(); ++i) {
      if(branch_start_nodes[i] == branch_end_nodes[i]) {
        std::ostringstream oss;
        oss<<"Found a branch with the same start and end node ("<<
          branch_start_nodes[i]<<"). Not allowed. ";
        throw std::logic_error(oss.str());
      }

      std::cout<<"Tree5"<<std::endl;

      auto it1 = map_node_to_id_.insert(
        std::pair<NodeType, uint>(branch_start_nodes[i], node_id_temp));

      if(it1.second) {
        // node encountered for the first time and inserted in the map_node_to_id_
        map_id_to_node_[node_id_temp] = branch_start_nodes[i];
        if(node_types[node_id_temp] == TIP) {
          node_types[node_id_temp] = INTERNAL;
        }
        branch_starts_temp[i] = node_id_temp;
        it_map_node_to_id_.push_back(it1.first);
        node_id_temp++;
      } else {
        // node encountered in a previous branch
        if(node_types[it1.first->second] == TIP) {
          // the previous encounter of the node was as a branch-end
          node_types[it1.first->second] = INTERNAL;
        } else {
          // do nothing
        }
        branch_starts_temp[i] = it1.first->second;
      }

      std::cout<<"Tree6"<<std::endl;

      auto it2 = map_node_to_id_.insert(std::pair<NodeType, uint>(branch_end_nodes[i], node_id_temp));

      std::cout<<"Tree6.1"<<std::endl;
      if(it2.second) {
        std::cout<<map_id_to_node_.size()<<" "<<node_id_temp<<" "
        <<branch_end_nodes.size()<<" "<<i<<std::endl;
        // node encountered for the first time and inserted in the map_node_to_id
        map_id_to_node_[node_id_temp] = branch_end_nodes[i];

        std::cout<<"Tree6.2"<<std::endl;
        if(node_types[node_id_temp] == ROOT) {
          // not known if the node has descendants, so we set its type to TIP.
          node_types[node_id_temp] = TIP;
        }
        std::cout<<"Tree6.3"<<std::endl;
        branch_ends_temp[i] = node_id_temp;
        ending_at[node_id_temp] = i;
        it_map_node_to_id_.push_back(it2.first);
        node_id_temp++;
        std::cout<<"Tree6.4"<<std::endl;
      } else {
        std::cout<<"Tree6.5"<<std::endl;
        // node has been previously encountered
        if(ending_at[it2.first->second] != NA_UINT) {
          std::cout<<"Tree6.6"<<std::endl;
          std::ostringstream oss;
          oss<<"Found at least two branches ending at the same node ("<<
            it2.first->first<<"). Check for cycles or repeated branches. ";
          throw std::logic_error(oss.str());
        } else {
          std::cout<<"Tree6.7"<<std::endl;
          if(node_types[it2.first->second] == ROOT) {
            // the previous enounters of the node were as branch-start -> set
            // the node's type to INTERNAL, because we know for sure that it
            // has descendants.
            std::cout<<"Tree6.8"<<std::endl;
            node_types[it2.first->second] = INTERNAL;
          }
          std::cout<<"Tree6.9"<<std::endl;
          branch_ends_temp[i] = it2.first->second;
          std::cout<<"Tree6.10"<<std::endl;
          ending_at[it2.first->second] = i;
          std::cout<<"Tree6.11"<<std::endl;
        }
      }
      std::cout<<"Tree7"<<std::endl;
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

    std::cout<<"Tree8"<<std::endl;
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
    std::cout<<"Tree9"<<std::endl;
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
      //map_node_to_id_[map_id_to_node_[i]] = node_ids[i];
      it_map_node_to_id_[i]->second = node_ids[i];
    }

    std::cout<<"Tree10"<<std::endl;

    this->map_id_to_node_ = At(map_id_to_node_, SortIndices(node_ids));

    this->id_parent_ = uvec(num_nodes_ - 1);

    if(branch_lengths.size() == num_nodes_ - 1) {
      this->lengths_ = std::vector<LengthType>(num_nodes_ - 1);
    } else if(branch_lengths.size() != 0) {
      std::ostringstream oss;
      oss<<"branch_lengths should be either empty or of size num_nodes_-1 ("<<
        num_nodes_-1<<") but is "<<branch_lengths.size()<<"."<<std::endl;
      throw std::invalid_argument(oss.str());
    }

    std::cout<<"Tree11"<<std::endl;

    if(HasBranchLengths()) {
      for(uint i = 0; i < num_nodes_ - 1; i++) {
        uint branch_start_i = node_ids[branch_starts_temp[i]];
        uint branch_end_i = node_ids[branch_ends_temp[i]];
        id_parent_[branch_end_i] = branch_start_i;
        lengths_[branch_end_i] = branch_lengths[i];
      }
      std::cout<<"Tree12"<<std::endl;
    } else {
      for(uint i = 0; i < num_nodes_ - 1; i++) {
        uint branch_start_i = node_ids[branch_starts_temp[i]];
        uint branch_end_i = node_ids[branch_ends_temp[i]];
        id_parent_[branch_end_i] = branch_start_i;
      }

      std::cout<<"Tree13"<<std::endl;
    }

    std::cout<<"Tree14"<<std::endl;
    init_id_child_nodes();
    std::cout<<"Tree15"<<std::endl;
  }

  uint num_nodes() const {
    return num_nodes_;
  }

  uint num_tips() const {
    return num_tips_;
  }

  bool HasBranchLengths() const {
    return lengths_.size() == id_parent_.size();
  }

  LengthType const& LengthOfBranch(uint i) const {
    if(i >= lengths_.size()) {
      std::ostringstream oss;
      oss<<"i is beyond the size of the lengths_ vector."<<
        "Check i and that the tree has branches."<<std::endl;
    }
    return lengths_[i];
  }

  std::vector<LengthType> const& lengths() const {
    return lengths_;
  }

  void SetLengthOfBranch(uint i, Length const& value) {
    if(!HasBranchLengths()) {
      std::ostringstream oss;
      oss<<"Trying to set a branch length on a tree without branch lengths. "<<
        "Use a SetBranchLengths method to add branch lengths first."<<std::endl;
      throw std::logic_error(oss.str());
    } else if(i >= lengths_.size()) {
      std::ostringstream oss;
      oss<<"i shoulc be smaller than "<<lengths_.size()<<" but was "<<i<<std::endl;
      throw std::out_of_range(oss.str());
    } else {
      lengths_[i] = value;
    }
  }

  void SetBranchLengths(std::vector<LengthType> const& lengths) {
    if(lengths.size() != 0 && lengths.size() != num_nodes_ - 1) {
      std::ostringstream oss;
      oss<<"lengths should be either empty or of size num_nodes_-1 ("<<
        num_nodes_-1<<") but is "<<lengths.size()<<"."<<std::endl;
    } else {
      lengths_ = lengths;
    }
  }

  void SetBranchLengths(std::vector<NodeType> const& nodes_branch_ends,
                        std::vector<LengthType> const& lengths) {
    if(nodes_branch_ends.size() != lengths.size()) {
      throw std::invalid_argument("The vectors nodes_branch_ends and lengths should be the same size.");
    }
    if( !HasBranchLengths() ) {
      if(nodes_branch_ends.size() != num_nodes_ - 1) {
        std::ostringstream oss;
        oss<<"Trying to set branch lengths on a tree without such."<<
          " In this case, the vectors nodes_branch_ends and lengths should have an"<<
            "element for each branch but their size is "<<
              nodes_branch_ends.size()<<" (should be "<<num_nodes_ - 1<<")."<<std::endl;
        throw std::invalid_argument(oss.str());
      }
      lengths_ = vec(num_nodes_ - 1);
    }
    std::vector<bool> visited(num_nodes_ - 1, false);
    for(uint i = 0; i < nodes_branch_ends.size(); ++i) {
      uint id = FindIdOfNode(nodes_branch_ends[i]);
      if(i == NA_UINT || i == num_nodes_ - 1) {
        std::ostringstream oss;
        oss<<"No branch ends at node identified as "<<id<<
          ". Check that nodes_branch_ends correspond to tips or internal nodes (excluding the root)"<<std::endl;
        throw std::logic_error(oss.str());
      } else if(visited[id]) {
        std::ostringstream oss;
        oss<<"Trying to set the length of the same branch twice. Check nodes_branch_ends for duplicates."<<std::endl;
        throw std::logic_error(oss.str());
      }
      visited[id] = true;
      lengths_[id] = lengths[i];
    }

  }

  // return the node at position id;
  NodeType const& FindNodeWithId(uint id) const {
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

  uvec const& FindChildren(uint i) const {
    if(i < this->num_tips()) {
      return EMPTY_UVEC;
    } else if(i - this->num_tips() < id_child_nodes_.size()) {
      return id_child_nodes_[i - this->num_tips()];
    } else {
      throw std::invalid_argument("i must be smaller than the number of nodes.");
    }
  }

  // returns a vector of positions in nodes in the order of their internally stored ids.
  uvec OrderNodes(std::vector<NodeType> const& nodes) const {
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


template<class Node, class Length>
class PruningTree: public Tree<Node, Length> {
public:
  typedef Node NodeType;
  typedef Length LengthType;

private:
  // default constructor;
  PruningTree() {}

protected:
  uvec ranges_id_visit_;
  uvec ranges_id_prune_;

public:

  PruningTree(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths):
  Tree<NodeType, LengthType>(branch_start_nodes, branch_end_nodes, branch_lengths),
  ranges_id_visit_(1, 0),
  ranges_id_prune_(1, 0) {

    std::cout<<"PruningTree1"<<std::endl;

    // insert a fictive branch leading to the root of the tree.
    uvec branch_ends = Seq(0, this->num_nodes_ - 1);

    uvec num_children_remaining(this->num_nodes_, 0);
    for(uint i : this->id_parent_) num_children_remaining[i]++;

    // start by pruning the tips of the tree
    uvec tips_this_level = Seq(0, this->num_tips_ - 1);

    uvec default_pos_vector;
    default_pos_vector.reserve(2);
    std::vector<uvec> pos_of_parent(this->num_nodes_ - this->num_tips_,
                                    default_pos_vector);

    uvec order_branches;
    order_branches.reserve(this->num_nodes_);

    std::cout<<"PruningTree2"<<std::endl;

    while(tips_this_level[0] != this->num_nodes_ - 1) {
      // while the root has not become a tip itself
      ranges_id_visit_.push_back(
        ranges_id_visit_[ranges_id_visit_.size() - 1] + tips_this_level.size());

      std::cout<<"PruningTree3"<<std::endl;

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

      std::cout<<"PruningTree4"<<std::endl;
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

        std::cout<<"PruningTree5"<<std::endl;
        ranges_id_prune_.push_back(
          ranges_id_prune_[ranges_id_prune_.size() - 1] + num_parent_updates);
      }

      std::cout<<"PruningTree6"<<std::endl;
      tips_this_level = tips_next_level;
    }

    if(this->HasBranchLengths()) {
      this->lengths_ = At(this->lengths_, order_branches);
    }

    std::cout<<"PruningTree7"<<std::endl;

    uvec id_old = order_branches;
    id_old.push_back(this->num_nodes_ - 1);

    this->id_parent_ = Match(At(this->id_parent_, order_branches),
                             id_old);

    std::cout<<"PruningTree8"<<std::endl;

    // update maps
    std::vector<NodeType> map_id_to_node(this->num_nodes_);
    for (uint i = 0; i < this->num_nodes_; i++) {
      map_id_to_node[i] = this->map_id_to_node_[id_old[i]];
      this->map_node_to_id_[map_id_to_node[i]] = i;
    }

    std::cout<<"PruningTree9"<<std::endl;

    std::swap(this->map_id_to_node_, map_id_to_node);

    std::cout<<"PruningTree10"<<std::endl;

    this->init_id_child_nodes();
    std::cout<<"PruningTree11"<<std::endl;
  }

  uint num_levels() const {
    return ranges_id_visit_.size() - 1;
  }

  uint num_parallel_ranges_prune() const {
    return ranges_id_prune_.size() - 1;
  }

  uvec const& ranges_id_visit() const {
    return ranges_id_visit_;
  }

  std::pair<uint, uint> RangeIdVisitNode(uint i_level) const {
    return std::pair<uint, uint>(ranges_id_visit_[i_level],
                                 ranges_id_visit_[i_level+1] - 1);
  }

  uvec const& ranges_id_prune() const {
    return ranges_id_prune_;
  }

  std::pair<uint, uint> RangeIdPruneNode(uint i_step) const {
    return std::pair<uint, uint>(ranges_id_prune_[i_step],
                                 ranges_id_prune_[i_step+1] - 1);
  }

  // A root-to-node distance vector in the order of pruning processing
  std::vector<LengthType> CalculateHeights(Length const& zero) const {
    std::vector<LengthType> h(this->num_nodes_, zero);
    for(int i = this->num_nodes_ - 2; i >= 0; i--) {
      h[i] = h[this->id_parent_[i]] + this->lengths_[i];
    }
    return h;
  }
};

enum ParallelMode {
  AUTO = 0,
  SINGLE_THREAD_LOOP_POSTORDER = 10,
  SINGLE_THREAD_LOOP_PRUNES = 11,
  SINGLE_THREAD_LOOP_VISITS = 12,
  MULTI_THREAD_LOOP_PRUNES = 21,
  MULTI_THREAD_LOOP_VISITS = 22,
  MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES = 23,
  MULTI_THREAD_VISIT_QUEUE = 24,
  HYBRID_LOOP_PRUNES = 31,
  HYBRID_LOOP_VISITS = 32,
  HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES = 33
};

inline std::ostream& operator<< (std::ostream& os, ParallelMode mode) {
  switch(mode) {
  case ParallelMode::AUTO: os<<"AUTO"; break;
  case ParallelMode::SINGLE_THREAD_LOOP_POSTORDER: os<<"SINGLE_THREAD_LOOP_POSTORDER"; break;
  case ParallelMode::SINGLE_THREAD_LOOP_PRUNES: os<<"SINGLE_THREAD_LOOP_PRUNES"; break;
  case ParallelMode::SINGLE_THREAD_LOOP_VISITS: os<<"SINGLE_THREAD_LOOP_VISITS"; break;
  case ParallelMode::MULTI_THREAD_LOOP_PRUNES: os<<"MULTI_THREAD_LOOP_PRUNES"; break;
  case ParallelMode::MULTI_THREAD_LOOP_VISITS: os<<"MULTI_THREAD_LOOP_VISITS"; break;
  case ParallelMode::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES: os<<"MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES"; break;
  case ParallelMode::MULTI_THREAD_VISIT_QUEUE: os<<"MULTI_THREAD_VISIT_QUEUE"; break;
  case ParallelMode::HYBRID_LOOP_PRUNES: os<<"HYBRID_LOOP_PRUNES"; break;
  case ParallelMode::HYBRID_LOOP_VISITS: os<<"HYBRID_LOOP_VISITS"; break;
  case ParallelMode::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES: os<<"HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES"; break;
  };
  return os<< static_cast<int>(mode);
}

template<class PruningSpec>
class ParallelPruning {
public:
  typedef ParallelMode PruningModeType;

protected:
  typedef typename PruningSpec::TreeType TreeType;

  TreeType const& ref_tree_;
  PruningSpec& ref_spec_;

  uint num_threads_;

  uint current_step_tuning_ = 0;
  uint fastest_step_tuning_ = 0;

  double min_duration_tuning_ = std::numeric_limits<double>::max();
  std::vector<double> durations_tuning_;

  const uvec min_sizes_chunk_ = {8}; //, 4, 8, 16, 32};

  const std::vector<ParallelMode> choices_mode_auto_ = {
    ParallelMode::SINGLE_THREAD_LOOP_POSTORDER,
    ParallelMode::SINGLE_THREAD_LOOP_PRUNES,
    ParallelMode::SINGLE_THREAD_LOOP_VISITS,
    ParallelMode::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES,
    ParallelMode::MULTI_THREAD_LOOP_VISITS,
    ParallelMode::MULTI_THREAD_VISIT_QUEUE
  };

  const std::vector<ParallelMode> choices_hybrid_mode_auto_ = {
    ParallelMode::HYBRID_LOOP_PRUNES,
    ParallelMode::HYBRID_LOOP_VISITS,
    ParallelMode::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES
  };

  class VisitQueue {
    std::mutex mutex_;
    std::condition_variable has_a_new_node_;

    TreeType const& ref_tree_;
    uvec queue_;
    uvec::iterator it_queue_begin;
    uvec::iterator it_queue_end;
    uvec num_non_visited_children_;
  public:

    // non-thread safe (called in single-thread mode)
    void Init(uvec const& num_children) {
      std::copy(num_children.begin(), num_children.end(),
                num_non_visited_children_.begin());
      it_queue_begin = queue_.begin();
      it_queue_end = queue_.begin() + ref_tree_.num_tips();
      std::iota(it_queue_begin, it_queue_end, 0);
    }

    bool IsTemporarilyEmpty() const {
      return it_queue_begin == it_queue_end & it_queue_end < queue_.end();
    }

    // thread-safe
    uint NextInQueue() {
      std::unique_lock<std::mutex> lock(mutex_);

      while( IsTemporarilyEmpty() ) {
        has_a_new_node_.wait(lock);
      }

      if(it_queue_begin < it_queue_end) {
        uint res = *it_queue_begin;
        ++it_queue_begin;
        return res;
      } else if(it_queue_begin == queue_.end()) {
        // algorithm thread should stop here. all waiting threads should be notified,
        // since no other elements will be inserted in the queue.
        has_a_new_node_.notify_all();
        return ref_tree_.num_nodes();
      } else {
        // algorithm thread continues to check for new node to visit
        // should never execute this
        //std::cout<<"Error returning NA_UINT from VisitQueue."<<std::endl;
        return NA_UINT;
      }
    }

    // thread-safe
    // if the parent of i becomes visit-able, it gets inserted in the
    // queue.
    void RemoveVisitedNode(uint i) {
      std::unique_lock<std::mutex> lock(mutex_);

      uint i_parent = ref_tree_.FindIdOfParent(i);
      num_non_visited_children_[i_parent - ref_tree_.num_tips()]--;
      if(num_non_visited_children_[i_parent - ref_tree_.num_tips()] == 0) {
        *it_queue_end = i_parent;
        *it_queue_end++;
        has_a_new_node_.notify_one();
      }
    }

    // non-thread-safe. should call Init() before using.
    VisitQueue(TreeType const& tree):
    ref_tree_(tree),
    queue_(tree.num_nodes()),
    it_queue_begin(queue_.begin()),
    it_queue_end(queue_.begin()),
    num_non_visited_children_(tree.num_nodes() - tree.num_tips()) {}

    // Copy initialization (non-thread-safe)
    VisitQueue(const VisitQueue& other): ref_tree_(other.ref_tree_) {
      auto other_begin = other.queue_.begin();
      queue_ = other.queue_;
      it_queue_begin = queue_.begin() + (other.it_queue_begin  - other_begin);
      it_queue_end = queue_.begin() + (other.it_queue_end  - other_begin);
      num_non_visited_children_ = other.num_non_visited_children_;
    }
  };

  uvec num_children_;
  VisitQueue visit_queue_;
public:

  ParallelPruning(TreeType const& tree, PruningSpec& spec):
  ref_tree_(tree),
  ref_spec_(spec),
  num_children_(tree.num_nodes() - tree.num_tips()),
  visit_queue_(tree) {

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

    for(uint i = tree.num_tips(); i < tree.num_nodes(); i++) {
      num_children_[i - tree.num_tips()] = tree.FindChildren(i).size();
    }
  }

  uint num_threads() const {
    return num_threads_;
  }

  bool IsTuning() const {
    return current_step_tuning_ < choices_mode_auto_.size() +
      min_sizes_chunk_.size() * choices_hybrid_mode_auto_.size();
  }

  uint VersionOPENMP() const {
#ifdef _OPENMP
    return _OPENMP;
#else
    return 0;
#endif
  }

  std::string ModeAutoCurrent() const {
    std::ostringstream oss;
    oss<<ModeAuto();
    return oss.str();
  }

  std::string ModeAutoStep(uint step) const {
    std::ostringstream oss;
    oss<<ModeAuto(step);
    return oss.str();
  }

  PruningModeType ModeAuto() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ModeAuto(step);
  }

  PruningModeType ModeAuto(uint step) const {
    if( step < choices_mode_auto_.size() ) {
      return choices_mode_auto_[step];
    } else {
      uint k = choices_hybrid_mode_auto_.size();
      uint l = step - choices_mode_auto_.size();
      return choices_hybrid_mode_auto_[(l/k) % k];
    }

  }

  uint min_size_chunk_visit() const {
    return min_sizes_chunk_[IndexMinSizeChunkVisit()];
  }

  uint min_size_chunk_prune() const {
    return min_sizes_chunk_[IndexMinSizeChunkPrune()];
  }

  uint fastest_step_tuning() const {
    return fastest_step_tuning_;
  }

  std::vector<double>  durations_tuning() const {
    return durations_tuning_;
  }

  void DoPruning(ParallelMode mode) {
    switch(mode) {
    case ParallelMode::SINGLE_THREAD_LOOP_POSTORDER: DoPruningSingleThreadLoopPostorder(); break;
    case ParallelMode::SINGLE_THREAD_LOOP_PRUNES: DoPruningSingleThreadLoopPrunes(); break;
    case ParallelMode::SINGLE_THREAD_LOOP_VISITS: DoPruningSingleThreadLoopVisits(); break;
    case ParallelMode::MULTI_THREAD_LOOP_PRUNES: DoPruningMultiThreadLoopPrunes(); break;
    case ParallelMode::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES: DoPruningMultiThreadLoopVisitsThenLoopPrunes(); break;
    case ParallelMode::MULTI_THREAD_LOOP_VISITS: DoPruningMultiThreadLoopVisits(); break;
    case ParallelMode::MULTI_THREAD_VISIT_QUEUE: DoPruningMultiThreadVisitQueue(); break;
    case ParallelMode::HYBRID_LOOP_PRUNES: DoPruningHybridLoopPrunes(); break;
    case ParallelMode::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES: DoPruningHybridLoopVisitsThenLoopPrunes(); break;
    case ParallelMode::HYBRID_LOOP_VISITS: DoPruningHybridLoopVisits(); break;
    default: DoPruningAuto();
    }
  }
protected:

  uint IndexMinSizeChunkVisit() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    //return (step / min_sizes_chunk_.size()) % min_sizes_chunk_.size();
    return step % min_sizes_chunk_.size();
  }

  uint IndexMinSizeChunkPrune() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return step % min_sizes_chunk_.size();
  }

  void DoPruningAuto() {

    std::chrono::steady_clock::time_point start, end;
    double duration;

    ParallelMode mode = ModeAuto();

    if( IsTuning() ) {

      start = std::chrono::steady_clock::now();
      DoPruning(mode);
      end = std::chrono::steady_clock::now();

      duration = std::chrono::duration<double, std::milli>(end - start).count();
      durations_tuning_.push_back(duration);
      if(duration < min_duration_tuning_) {
        min_duration_tuning_ = duration;
        fastest_step_tuning_ = current_step_tuning_;
      }
      current_step_tuning_++;

    } else {
      DoPruning(mode);
    }

  }

  void DoPruningSingleThreadLoopPostorder() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }

    for(uint i = 0; i < ref_tree_.num_nodes() - 1; i++) {
      ref_spec_.VisitNode(i);
      ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
    }
  }

  void DoPruningSingleThreadLoopPrunes() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }

    for(int i_prune = 0; i_prune < ref_tree_.num_parallel_ranges_prune(); i_prune++) {
      auto range_prune = ref_tree_.RangeIdPruneNode(i_prune);

    _PRAGMA_OMP_SIMD
      for(uint i = range_prune.first; i <= range_prune.second; i++) {
        ref_spec_.VisitNode(i);
        ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
      }
    }
  }

  void DoPruningSingleThreadLoopVisits() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }

    for(int i_level = 0; i_level < ref_tree_.num_levels(); i_level++) {
      auto range_visit = ref_tree_.RangeIdVisitNode(i_level);
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit.first; i <= range_visit.second; i++) {
        if(i < ref_tree_.num_tips()) {
          // i is a tip (only Visit)
          ref_spec_.VisitNode(i);
        } else {
          // i is internal or root
          for(uint j: ref_tree_.FindChildren(i)) {
            ref_spec_.PruneNode(j, i);
          }
          ref_spec_.VisitNode(i);
        }
      }
    }

    // VisitNode not called on the root node
    for(uint j: ref_tree_.FindChildren(ref_tree_.num_nodes() - 1)) {
      ref_spec_.PruneNode(j, ref_tree_.num_nodes() - 1);
    }
  }

  void DoPruningMultiThreadLoopVisitsThenLoopPrunes() {

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
    ref_spec_.InitNode(i);
  }

  uint i_prune = 0;
  for(uint i_level = 0; i_level < ref_tree_.num_levels(); i_level++) {

#pragma omp barrier

    auto range_visit = ref_tree_.RangeIdVisitNode(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit.first; i <= range_visit.second; i++) {
        ref_spec_.VisitNode(i);
      }

      uint num_branches_done = 0;

    while(num_branches_done != range_visit.second - range_visit.first + 1) {
#pragma omp barrier
      auto range_prune = ref_tree_.RangeIdPruneNode(i_prune);

      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_prune.first; i <= range_prune.second; i++) {
          ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
        }

        num_branches_done +=  range_prune.second - range_prune.first + 1;
      ++i_prune;
    }
  }
}
  }

  void DoPruningMultiThreadLoopVisits() {
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }

    for(int i_level = 0; i_level < ref_tree_.num_levels(); i_level++) {
      auto range_visit = ref_tree_.RangeIdVisitNode(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit.first; i <= range_visit.second; i++) {
        if(i < ref_tree_.num_tips()) {
          // i is a tip (only Visit)
          ref_spec_.VisitNode(i);
        } else {
          // i is internal or root
          for(uint j: ref_tree_.FindChildren(i)) {
            ref_spec_.PruneNode(j, i);
          }
          ref_spec_.VisitNode(i);
        }
      }
    }
}
    // VisitNode not called on the root node
    for(uint j: ref_tree_.FindChildren(ref_tree_.num_nodes() - 1)) {
      ref_spec_.PruneNode(j, ref_tree_.num_nodes() - 1);
    }
  }

  void DoPruningMultiThreadVisitQueue() {
    visit_queue_.Init(num_children_);
#pragma omp parallel
{
  while(true) {
    uint i = visit_queue_.NextInQueue();
    if(i == NA_UINT) {
      continue;
    } else if(i == ref_tree_.num_nodes()) {
      break;
    } else if(i < ref_tree_.num_tips()) {
      // i is a tip (only Visit)
      ref_spec_.InitNode(i);
      ref_spec_.VisitNode(i);
      visit_queue_.RemoveVisitedNode(i);
    } else if(i < ref_tree_.num_nodes() - 1){
      // i is internal
      ref_spec_.InitNode(i);
      uvec const& children = ref_tree_.FindChildren(i);
      for(uint j: children) {
        ref_spec_.PruneNode(j, i);
      }
      ref_spec_.VisitNode(i);
      visit_queue_.RemoveVisitedNode(i);
    } else {
      // i is the root
      ref_spec_.InitNode(i);
      uvec const& children = ref_tree_.FindChildren(i);
      for(uint j: children) {
        ref_spec_.PruneNode(j, i);
      }
      // don't visit the root
    }
  }
}
  }

  void DoPruningMultiThreadLoopPrunes() {

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
    ref_spec_.InitNode(i);
  }

  for(int i_prune = 0; i_prune < ref_tree_.num_parallel_ranges_prune(); i_prune++) {
    auto range_prune = ref_tree_.RangeIdPruneNode(i_prune);

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_prune.first; i <= range_prune.second; i++) {
        ref_spec_.VisitNode(i);
        ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
      }
  }
}
  }

  void DoPruningHybridLoopVisitsThenLoopPrunes() {
    uint min_size_chunk_visit = this->min_size_chunk_visit();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }

    uint i_prune = 0;
  for(int i_level = 0; i_level < ref_tree_.num_levels(); i_level++) {
    auto range_visit = ref_tree_.RangeIdVisitNode(i_level);
#pragma omp barrier
    if(range_visit.second - range_visit.first + 1 >
        num_threads_ * min_size_chunk_visit) {
      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_visit.first; i <= range_visit.second; i++) {
          ref_spec_.VisitNode(i);
        }
    } else if(tid == 0) {
      // only the master thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit.first; i <= range_visit.second; i++) {
        ref_spec_.VisitNode(i);
      }
    }

    if (tid == 0) {
      // only one (master) thread executes this
      uint num_branches_done = 0;
      while(num_branches_done != range_visit.second - range_visit.first + 1) {
        auto range_prune = ref_tree_.RangeIdPruneNode(i_prune);
        _PRAGMA_OMP_SIMD
          for(uint i = range_prune.first; i <= range_prune.second; i++) {
            ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
          }

          num_branches_done +=  range_prune.second - range_prune.first + 1;
        ++i_prune;
      }
    }
  }
}
  }

  void DoPruningHybridLoopPrunes() {
    uint min_size_chunk_prune = this->min_size_chunk_prune();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }


  for(int i_prune = 0; i_prune < ref_tree_.num_parallel_ranges_prune(); i_prune++) {
      auto range_prune = ref_tree_.RangeIdPruneNode(i_prune);
#pragma omp barrier
      if (range_prune.second - range_prune.first + 1 >
            num_threads_ * min_size_chunk_prune) {
        _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_prune.first; i <= range_prune.second; i++) {
          ref_spec_.VisitNode(i);
          ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
        }
      } else if (tid == 0) {
        // only one (master) thread executes this
        _PRAGMA_OMP_SIMD
        for(uint i = range_prune.first; i <= range_prune.second; i++) {
          ref_spec_.VisitNode(i);
          ref_spec_.PruneNode(i, ref_tree_.FindIdOfParent(i));
        }
      }
    }
}
  }

  void DoPruningHybridLoopVisits() {
    uint min_size_chunk_visit = this->min_size_chunk_visit();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ref_tree_.num_nodes(); i++) {
      ref_spec_.InitNode(i);
    }

  for(int i_level = 0; i_level < ref_tree_.num_levels(); i_level++) {
    auto range_visit = ref_tree_.RangeIdVisitNode(i_level);
#pragma omp barrier
    if(range_visit.second - range_visit.first + 1 >
         num_threads_ * min_size_chunk_visit) {
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit.first; i <= range_visit.second; i++) {
        if(i < ref_tree_.num_tips()) {
          // i is a tip (only Visit)
          ref_spec_.VisitNode(i);
        } else if(i < ref_tree_.num_nodes() - 1){
          // i is internal
          for(uint j: ref_tree_.FindChildren(i)) {
            ref_spec_.PruneNode(j, i);
          }
          ref_spec_.VisitNode(i);
        }
      }
    } else if(tid == 0) {
      // only the master thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit.first; i <= range_visit.second; i++) {
        if(i < ref_tree_.num_tips()) {
          // i is a tip (only Visit)
          ref_spec_.VisitNode(i);
        } else if(i < ref_tree_.num_nodes() - 1){
          // i is internal
          for(uint j: ref_tree_.FindChildren(i)) {
            ref_spec_.PruneNode(j, i);
          }
          ref_spec_.VisitNode(i);
        }
      }
    }
  }
}
    // VisitNode not called on the root
    for(uint j: ref_tree_.FindChildren(ref_tree_.num_nodes() - 1)) {
      ref_spec_.PruneNode(j, ref_tree_.num_nodes() - 1);
    }
  }

};


template<class PruningSpec>
class PruningTask {
public:
  typedef PruningSpec PruningSpecType;
  typedef typename PruningSpec::TreeType TreeType;
  typedef typename PruningSpec::PruningAlgorithmType PruningAlgorithmType;
  typedef typename PruningAlgorithmType::PruningModeType PruningModeType;
  typedef typename TreeType::NodeType NodeType;
  typedef typename TreeType::LengthType LengthType;
  typedef typename PruningSpecType::InputDataType InputDataType;
  typedef typename PruningSpecType::ParameterType ParameterType;
  typedef typename PruningSpecType::NodeStateType NodeStateType;

  PruningTask(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths,
    InputDataType const& data):
  tree_(branch_start_nodes, branch_end_nodes, branch_lengths),
  spec_(tree_, data),
  algorithm_(tree_, spec_) {}

  NodeStateType DoPruning(ParameterType const& par, uint mode) {
    spec_.SetParameter(par);
    algorithm_.DoPruning(static_cast<PruningModeType>(mode));
    return spec_.StateAtRoot();
  }

  TreeType & tree() {
    return tree_;
  }
  PruningSpec & spec() {
    return spec_;
  }
  PruningAlgorithmType & algorithm() {
    return algorithm_;
  }
private:
  TreeType tree_;
  PruningSpec spec_;
  PruningAlgorithmType algorithm_;
};


// implicit interface
template<class Tree> class PruningSpecification {
protected:
  Tree const& ref_tree_;
  PruningSpecification(Tree const& tree): ref_tree_(tree) {}
};
}
#endif // ParallelPruning_ParallelPruning_H_
