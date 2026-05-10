//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// Bridge between the construction-time BuilderNode tree (tree_builder.h)
// and the query-time compact PDLTreeCore (tree_core.h). This is the only
// header that pulls in both: tree_builder.h stays SDSL-free for fast unit
// tests, and tree_core.h stays builder-free for disk-loaded indexes.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

#include "storage_policy.h"
#include "tree_builder.h"
#include "tree_core.h"

namespace dret::pdl {

// Walk a post-Task-11 BuilderNode tree (every node has node_id assigned in
// post-order) and assemble the compact PDLTreeCore. The sentinel value
// n_nodes is used in first_child / next_sibling vectors to mark "none";
// bit_compress on the populated int_vector picks the smallest width that
// still fits the sentinel.
template <typename TBitvector, typename TBvRank, typename TBvSelect,
          typename TIntVector, typename TStoredSetCodec>
void BuildPDLTreeCoreFromBuilder(
    PDLTreeCore<TBitvector, TBvRank, TBvSelect, TIntVector, TStoredSetCodec>& t_core,
    const BuilderNode* t_root,
    std::size_t t_n_nodes,
    std::size_t t_n_doc,
    uint32_t t_block_size,
    float t_storing_factor,
    StoragePolicy t_policy) {
  if (!t_root || t_n_nodes == 0) {
    t_core.Assemble(TIntVector{}, TIntVector{}, TIntVector{}, TIntVector{},
                    TBitvector{}, TStoredSetCodec{},
                    t_n_doc, t_block_size, t_storing_factor, t_policy);
    return;
  }

  TIntVector node_starts(t_n_nodes, 0);
  TIntVector node_ends(t_n_nodes, 0);
  TIntVector first_child(t_n_nodes, t_n_nodes);
  TIntVector next_sibling(t_n_nodes, t_n_nodes);
  sdsl::bit_vector selected_bv(t_n_nodes, 0);

  std::vector<const BuilderNode*> stack{t_root};
  while (!stack.empty()) {
    const auto* n = stack.back();
    stack.pop_back();
    const std::size_t id = n->node_id;
    node_starts[id] = n->sp;
    node_ends[id] = n->ep;
    first_child[id] = n->first_child ? n->first_child->node_id : t_n_nodes;
    next_sibling[id] = n->next_sibling ? n->next_sibling->node_id : t_n_nodes;
    if (n->selected) selected_bv[id] = 1;
    for (const auto* c = n->first_child; c; c = c->next_sibling) {
      stack.push_back(c);
    }
  }

  sdsl::util::bit_compress(node_starts);
  sdsl::util::bit_compress(node_ends);
  sdsl::util::bit_compress(first_child);
  sdsl::util::bit_compress(next_sibling);

  TBitvector selected_marker(selected_bv);

  t_core.Assemble(std::move(node_starts), std::move(node_ends),
                  std::move(first_child), std::move(next_sibling),
                  std::move(selected_marker), TStoredSetCodec{},
                  t_n_doc, t_block_size, t_storing_factor, t_policy);
}

}  // namespace dret::pdl
