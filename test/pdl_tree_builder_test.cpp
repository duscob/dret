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
using dret::pdl::CollapseSubtreesByBlockSize;
using dret::pdl::InsertExplicitLeaves;

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

// Block-size collapse on the same lcp = [0, 1, 4, 2, 1, 5, 3, 0] (n=7) tree
// used by CollapsePruning: tree builds as
//   root([0,7))
//     A([0,7), d=1)
//       C([1,4), d=2)
//         B([1,3), d=4)
//       E([4,7), d=3)
//         D([4,6), d=5)
// We trust that shape across the collapse tests below. Re-trace the LCP
// stack algorithm on paper before changing anything.
auto BuildCollapseFixture(BuilderPool& pool) {
  std::vector<std::size_t> lcp = {0, 1, 4, 2, 1, 5, 3, 0};
  return BuildSparseSuffixTree(pool, LcpFromVector(lcp), 7);
}

TEST(PDLTreeBuilder, CollapseNullRootIsNoop) {
  CollapseSubtreesByBlockSize(nullptr, 100);
}

TEST(PDLTreeBuilder, CollapseZeroBlockSizeChangesNothing) {
  BuilderPool pool;
  auto* root = BuildCollapseFixture(pool);

  CollapseSubtreesByBlockSize(root, 0);

  // Whole tree unchanged: A still has C and E; C still has B; E still has D.
  ASSERT_NE(root->first_child, nullptr);
  auto* a = root->first_child;
  ASSERT_NE(a->first_child, nullptr);
  auto* c = a->first_child;
  auto* e = c->next_sibling;
  ASSERT_NE(e, nullptr);
  EXPECT_NE(c->first_child, nullptr);
  EXPECT_NE(e->first_child, nullptr);
}

TEST(PDLTreeBuilder, CollapsePrunesGrandchildrenAtThreshold) {
  BuilderPool pool;
  auto* root = BuildCollapseFixture(pool);

  // C and E each have length 3; collapse with block_size=3 prunes their
  // children (B and D) but keeps C and E as childless nodes.
  CollapseSubtreesByBlockSize(root, 3);

  auto* a = root->first_child;
  ASSERT_NE(a, nullptr);
  auto* c = a->first_child;
  ASSERT_NE(c, nullptr);
  auto* e = c->next_sibling;
  ASSERT_NE(e, nullptr);
  EXPECT_EQ(c->first_child, nullptr) << "C should be collapsed";
  EXPECT_EQ(e->first_child, nullptr) << "E should be collapsed";
  // Sibling link between C and E preserved; both still cover their original
  // SA ranges (boundary info needed for raw-range fallback).
  EXPECT_EQ(c->sp, 1u);
  EXPECT_EQ(c->ep, 4u);
  EXPECT_EQ(e->sp, 4u);
  EXPECT_EQ(e->ep, 7u);
}

TEST(PDLTreeBuilder, CollapseAtRootKillsWholeTree) {
  BuilderPool pool;
  auto* root = BuildCollapseFixture(pool);

  // Root has length 7; block_size=7 collapses at root.
  CollapseSubtreesByBlockSize(root, 7);

  EXPECT_EQ(root->first_child, nullptr);
  EXPECT_EQ(root->sp, 0u);
  EXPECT_EQ(root->ep, 7u);
}

TEST(PDLTreeBuilder, CollapseStopsAtFirstSmallAncestor) {
  BuilderPool pool;
  auto* root = BuildCollapseFixture(pool);

  // block_size=2: leaves B and D have length 2 and get touched (no children
  // to remove); C and E have length 3 (>2), so descend into them. Result:
  // tree shape unchanged.
  CollapseSubtreesByBlockSize(root, 2);

  auto* a = root->first_child;
  ASSERT_NE(a, nullptr);
  auto* c = a->first_child;
  auto* e = c->next_sibling;
  ASSERT_NE(c, nullptr);
  ASSERT_NE(e, nullptr);
  ASSERT_NE(c->first_child, nullptr);
  ASSERT_NE(e->first_child, nullptr);
  EXPECT_EQ(c->first_child->sp, 1u);
  EXPECT_EQ(c->first_child->ep, 3u);
  EXPECT_EQ(e->first_child->sp, 4u);
  EXPECT_EQ(e->first_child->ep, 6u);
}

// Walks the tree and collects [sp, ep) of every "coverable" node — i.e.
// every node that has no children (collapsed blocks AND explicit leaves).
// Returned in left-to-right SA order via DFS that visits children in
// sibling order.
std::vector<std::pair<std::size_t, std::size_t>> CoverableRanges(const BuilderNode* root) {
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  // Use a recursive lambda for natural DFS order.
  auto visit = [&](auto& self, const BuilderNode* n) -> void {
    if (!n->first_child) {
      ranges.emplace_back(n->sp, n->ep);
      return;
    }
    for (auto* c = n->first_child; c; c = c->next_sibling) self(self, c);
  };
  visit(visit, root);
  return ranges;
}

// Acceptance criterion (Task 8): every SA position belongs to a leaf-like
// coverable unit. Verified by asserting CoverableRanges partitions [0, n).
TEST(PDLTreeBuilder, InsertLeavesNullRootIsNoop) {
  BuilderPool pool;
  InsertExplicitLeaves(pool, nullptr);
}

TEST(PDLTreeBuilder, NoLeavesNeededWhenChildrenAlreadyPartition) {
  // "abab" tree: root has two children whose ranges already cover [0, 4).
  std::vector<std::size_t> lcp = {0, 2, 0, 1, 0};
  BuilderPool pool;
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 4);
  std::size_t before = pool.size();

  InsertExplicitLeaves(pool, root);

  EXPECT_EQ(pool.size(), before) << "no leaves should be created";
  EXPECT_THAT(CoverableRanges(root),
              testing::ElementsAre(std::make_pair(0u, 2u), std::make_pair(2u, 4u)));
}

TEST(PDLTreeBuilder, LeavesFillGapAtFrontOfChain) {
  // "aaa" tree: root -> A([0,3)) -> B([1,3)). A's children only cover
  // [1, 3), so position 0 needs an explicit leaf as A's first child.
  std::vector<std::size_t> lcp = {0, 1, 2, 0};
  BuilderPool pool;
  auto* root = BuildSparseSuffixTree(pool, LcpFromVector(lcp), 3);

  InsertExplicitLeaves(pool, root);

  // Coverable ranges should be {[0,1), [1,3)} — the explicit leaf and B.
  EXPECT_THAT(CoverableRanges(root),
              testing::ElementsAre(std::make_pair(0u, 1u), std::make_pair(1u, 3u)));

  // The leaf is is_explicit_leaf=true and has depth 0 (drl convention).
  auto* a = root->first_child;
  ASSERT_NE(a, nullptr);
  auto* leaf = a->first_child;
  ASSERT_NE(leaf, nullptr);
  EXPECT_TRUE(leaf->is_explicit_leaf);
  EXPECT_EQ(leaf->depth, 0u);
  EXPECT_EQ(leaf->sp, 0u);
  EXPECT_EQ(leaf->ep, 1u);
  // Sibling chain leaf -> B intact.
  ASSERT_NE(leaf->next_sibling, nullptr);
  EXPECT_EQ(leaf->next_sibling->sp, 1u);
  EXPECT_EQ(leaf->next_sibling->ep, 3u);
}

TEST(PDLTreeBuilder, LeavesFillTrailingAndInteriorGaps) {
  // The collapse fixture without collapse: tree has gaps under A
  // (position 0 before C), under C (position 3 between B and C's end),
  // and under E (position 6 between D and E's end).
  BuilderPool pool;
  auto* root = BuildCollapseFixture(pool);

  InsertExplicitLeaves(pool, root);

  // Coverable ranges in SA order partition [0, 7).
  auto ranges = CoverableRanges(root);
  ASSERT_FALSE(ranges.empty());
  EXPECT_EQ(ranges.front().first, 0u);
  EXPECT_EQ(ranges.back().second, 7u);
  for (std::size_t i = 1; i < ranges.size(); ++i) {
    EXPECT_EQ(ranges[i - 1].second, ranges[i].first)
        << "gap between coverable range " << i - 1 << " and " << i;
  }
}

TEST(PDLTreeBuilder, CollapsedBlockSkippedAndStillCoversItsRange) {
  // After collapse with block_size=3 the C and E nodes lose their
  // grandchildren. addLeaves must not insert leaves under them — the block
  // covers its own range via raw-fallback.
  BuilderPool pool;
  auto* root = BuildCollapseFixture(pool);
  CollapseSubtreesByBlockSize(root, 3);

  InsertExplicitLeaves(pool, root);

  // After: A has [explicit leaf at pos 0, C, E]. C and E have no children.
  // Coverable ranges partition [0, 7).
  auto ranges = CoverableRanges(root);
  EXPECT_THAT(ranges, testing::ElementsAre(std::make_pair(0u, 1u),
                                           std::make_pair(1u, 4u),
                                           std::make_pair(4u, 7u)));
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
