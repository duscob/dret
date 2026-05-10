//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// Throwaway construction data for the PDL sparse suffix tree. These structs
// only exist while the tree is being built (Tasks 6-11); after that the
// tree is compacted into PDLTreeCore's serialized form and the builder is
// discarded. They must not appear in query-time headers.
//
// The linked-list shape (first_child / next_sibling) mirrors drl's
// PDLTreeNode (drl/include/drl/pdltree.h:192-238) so the stack/LCP port in
// Task 6 stays close to the reference implementation. Differences from drl:
// SA ranges are half-open [sp, ep) (drl uses inclusive); the explicit
// `selected` flag captures the storage-policy decision from Task 10.
//

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace dret::pdl {

struct BuilderNode {
  static constexpr std::size_t kInvalidId = std::numeric_limits<std::size_t>::max();

  // Half-open SA interval [sp, ep). Leaves have ep == sp + 1.
  std::size_t sp = 0;
  std::size_t ep = 0;

  // String depth (LCP value at the open of this interval). Leaves keep the
  // depth of the deepest internal ancestor; the algorithm in Task 6 sets it
  // when the leaf is created.
  std::size_t depth = 0;

  // Tree links. Linked-list children (first_child + next_sibling) keep the
  // construction algorithm allocation-light; a parent owns its first child
  // and that child owns its sibling chain.
  BuilderNode* parent = nullptr;
  BuilderNode* first_child = nullptr;
  BuilderNode* next_sibling = nullptr;

  // Document set computed bottom-up in Task 9. Empty when contains_all is
  // true; Task 9 also collapses full sets into the all-doc sentinel.
  std::vector<std::size_t> docs;
  bool contains_all = false;

  // Set by the storage policy in Task 10. Selected nodes get a stable
  // node_id (Task 11) and end up in PDLTreeCore's stored-set codec.
  bool selected = false;

  // Explicit leaves added by Task 8 to cover SA positions not represented
  // by internal LCP nodes after collapse. Distinguishing leaves from
  // single-position internals matters for some storage policies and for
  // navigation later on.
  bool is_explicit_leaf = false;

  // Weighted size used by the OriginalDrl selection rule in Task 10:
  // sum of stored-set sizes from the subtree rooted here. 0 until Task 10.
  std::size_t stored_documents = 0;

  // Final id assigned in Task 11. kInvalidId until then.
  std::size_t node_id = kInvalidId;

  BuilderNode() = default;
  BuilderNode(std::size_t t_sp, std::size_t t_ep, std::size_t t_depth)
      : sp(t_sp), ep(t_ep), depth(t_depth) {}
};

// Arena ownership for builder nodes; pointers stay stable for the lifetime
// of the pool. Discard the pool after compacting into PDLTreeCore.
class BuilderPool {
 public:
  BuilderNode* create() {
    auto node = std::make_unique<BuilderNode>();
    auto* raw = node.get();
    nodes_.push_back(std::move(node));
    return raw;
  }

  BuilderNode* create(std::size_t t_sp, std::size_t t_ep, std::size_t t_depth) {
    auto node = std::make_unique<BuilderNode>(t_sp, t_ep, t_depth);
    auto* raw = node.get();
    nodes_.push_back(std::move(node));
    return raw;
  }

  std::size_t size() const { return nodes_.size(); }

 private:
  std::vector<std::unique_ptr<BuilderNode>> nodes_;
};

// Append child as the first child of parent (drl's addChild order). Cheap
// and order-preserving in reverse; the construction stack will produce
// children in left-to-right order, so callers wanting that order should
// reverse the chain after building or use appendChild below.
inline void prependChild(BuilderNode* parent, BuilderNode* child) {
  child->parent = parent;
  child->next_sibling = parent->first_child;
  parent->first_child = child;
}

// Append child as the last child of parent. O(siblings) but keeps the
// natural left-to-right order produced by the LCP stack algorithm.
inline void appendChild(BuilderNode* parent, BuilderNode* child) {
  child->parent = parent;
  child->next_sibling = nullptr;
  if (!parent->first_child) {
    parent->first_child = child;
    return;
  }
  auto* tail = parent->first_child;
  while (tail->next_sibling) tail = tail->next_sibling;
  tail->next_sibling = child;
}

// Build a sparse suffix tree from an LCP array via the classical stack
// algorithm. Ports drl's loop at drl/src/pdltree.cpp:336-358 and converts
// drl's inclusive [sp, ep] ranges to dret's half-open [sp, ep) at the
// boundary (drl: range.second = i - 1; dret: ep = i).
//
// `t_lcp` must be indexable for i in [0, t_n] with t_lcp(0) == 0 and
// t_lcp(t_n) == 0 (sentinel forcing the stack to flush). For i in
// [1, t_n - 1], t_lcp(i) == LCP(SA[i-1], SA[i]).
//
// Returns the root node (depth 0, covering [0, t_n)) allocated from
// t_pool. Internal nodes are added below the root. Block-size collapse
// (Task 7) and explicit-leaf insertion (Task 8) are NOT part of this
// function — they are applied after construction. Returns nullptr for
// t_n == 0.
template <typename TLCPFn>
BuilderNode* BuildSparseSuffixTree(BuilderPool& t_pool, TLCPFn&& t_lcp, std::size_t t_n) {
  if (t_n == 0) return nullptr;

  std::vector<BuilderNode*> stack;
  stack.reserve(64);
  stack.push_back(t_pool.create(0, 0, 0));
  BuilderNode* root = stack.back();
  BuilderNode* prev = nullptr;

  for (std::size_t i = 1; i <= t_n; ++i) {
    std::size_t left = i - 1;
    const std::size_t lcp_i = static_cast<std::size_t>(t_lcp(i));

    while (lcp_i < stack.back()->depth) {
      stack.back()->ep = i;
      prev = stack.back();
      stack.pop_back();
      root = prev;
      left = prev->sp;
      if (lcp_i <= stack.back()->depth) {
        appendChild(stack.back(), prev);
        prev = nullptr;
      }
    }

    if (lcp_i > stack.back()->depth) {
      auto* curr = t_pool.create(left, left, lcp_i);
      if (prev) {
        appendChild(curr, prev);
        prev = nullptr;
      }
      stack.push_back(curr);
    }
  }

  while (root->parent) root = root->parent;
  root->ep = t_n;
  return root;
}

}  // namespace dret::pdl
