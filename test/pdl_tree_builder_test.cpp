//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// Acceptance tests for Task 6 (LCP-based sparse suffix-tree construction)
// from docs/pdl_indexes_tasks.md. The hand-traced expected trees below
// match what the classical stack/LCP algorithm produces on the listed
// LCP arrays; trace them again before changing any expected value.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "dret/pdl/tree_builder.h"

namespace {

using dret::pdl::BuilderNode;
using dret::pdl::BuilderPool;
using dret::pdl::BuildSparseSuffixTree;

struct NodeView {
  std::size_t sp = 0;
  std::size_t ep = 0;
  std::size_t depth = 0;
  std::vector<NodeView> children;

  bool operator==(const NodeView& o) const {
    return sp == o.sp && ep == o.ep && depth == o.depth && children == o.children;
  }
};

NodeView Capture(const BuilderNode* n) {
  NodeView v{n->sp, n->ep, n->depth, {}};
  for (auto* c = n->first_child; c; c = c->next_sibling) {
    v.children.push_back(Capture(c));
  }
  return v;
}

auto LcpFromVector(const std::vector<std::size_t>& v) {
  return [&v](std::size_t i) { return v.at(i); };
}

// "aaa" — sorted suffixes "a","aa","aaa". lcp = [0, 1, 2, 0].
// Expected: a single chain root -> A(depth 1) -> B(depth 2).
TEST(PDLTreeBuilder, RepeatedAaa) {
  std::vector<std::size_t> lcp = {0, 1, 2, 0};
  BuilderPool pool;
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 3);
  ASSERT_NE(root, nullptr);

  NodeView expected{
      0, 3, 0,
      {{0, 3, 1, {{1, 3, 2, {}}}}},
  };
  EXPECT_EQ(Capture(root), expected);
}

// "abab" — sorted suffixes "ab","abab","b","bab". lcp = [0, 2, 0, 1, 0].
// Expected: root has two children — A(depth 2, [0,2)) and B(depth 1, [2,4)).
TEST(PDLTreeBuilder, BranchingAbab) {
  std::vector<std::size_t> lcp = {0, 2, 0, 1, 0};
  BuilderPool pool;
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 4);
  ASSERT_NE(root, nullptr);

  NodeView expected{
      0, 4, 0,
      {
          {0, 2, 2, {}},
          {2, 4, 1, {}},
      },
  };
  EXPECT_EQ(Capture(root), expected);
}

// Single suffix — lcp = [0, 0]. Root only.
TEST(PDLTreeBuilder, SinglePosition) {
  std::vector<std::size_t> lcp = {0, 0};
  BuilderPool pool;
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 1);
  ASSERT_NE(root, nullptr);

  NodeView expected{0, 1, 0, {}};
  EXPECT_EQ(Capture(root), expected);
}

// Empty input — returns null.
TEST(PDLTreeBuilder, Empty) {
  BuilderPool pool;
  std::vector<std::size_t> lcp = {0};
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 0);
  EXPECT_EQ(root, nullptr);
}

// Half-open intervals: every internal node's [sp, ep) lies inside its
// parent's range, and sibling ranges are disjoint and ordered.
TEST(PDLTreeBuilder, HalfOpenInvariants) {
  // "abcabc" — sorted suffixes "abc","abcabc","bc","bcabc","c","cabc".
  // lcp = [0, 3, 0, 2, 0, 1, 0].
  std::vector<std::size_t> lcp = {0, 3, 0, 2, 0, 1, 0};
  BuilderPool pool;
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 6);
  ASSERT_NE(root, nullptr);

  std::vector<const BuilderNode*> stack{root};
  while (!stack.empty()) {
    const auto* n = stack.back();
    stack.pop_back();
    ASSERT_LT(n->sp, n->ep) << "node has empty half-open range";
    const BuilderNode* prev = nullptr;
    for (const auto* c = n->first_child; c; c = c->next_sibling) {
      EXPECT_GE(c->sp, n->sp);
      EXPECT_LE(c->ep, n->ep);
      EXPECT_GT(c->depth, n->depth);
      if (prev) EXPECT_GE(c->sp, prev->ep) << "siblings overlap";
      prev = c;
      stack.push_back(c);
    }
  }
}

}  // namespace
