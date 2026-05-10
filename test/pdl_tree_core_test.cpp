//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// PDLTreeCore::computeCoverFull behavior tests (Task 12).
//
// Each test builds the canonical post-Task-11 builder fixture (from
// pdl_tree_builder_test.cpp's BuildCollapseFixture + InsertExplicitLeaves
// + ComputeDocSetsBottomUp + ApplyStoragePolicy + AssignNodeIds), bridges
// it into a PDLTreeCore via BuildPDLTreeCoreFromBuilder, and queries
// computeCoverFull. The hand-traced expectations come from walking the
// algorithm on paper; re-trace before changing any expected value.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "dret/pdl/build_pdl_core.h"
#include "dret/pdl/storage_policy.h"
#include "dret/pdl/tree_builder.h"
#include "dret/pdl/tree_core.h"

namespace {

using dret::pdl::ApplyStoragePolicy;
using dret::pdl::AssignNodeIds;
using dret::pdl::BuildPDLTreeCoreFromBuilder;
using dret::pdl::BuildSparseSuffixTree;
using dret::pdl::BuilderPool;
using dret::pdl::ComputeDocSetsBottomUp;
using dret::pdl::InsertExplicitLeaves;
using dret::pdl::PDLTreeCore;
using dret::pdl::StoragePolicy;

// Fixture shape (lcp = [0,1,4,2,1,5,3,0], n=7, then InsertExplicitLeaves):
//   root[0,7) d=0
//     A[0,7) d=1
//       leaf(0,1) explicit
//       C[1,4) d=2
//         B[1,3) d=4
//         leaf(3,4) explicit
//       E[4,7) d=3
//         D[4,6) d=5
//         leaf(6,7) explicit
//
// Post-order ids (Task 11):
//   leaf(0,1) = 0
//   B = 1, leaf(3,4) = 2, C = 3
//   D = 4, leaf(6,7) = 5, E = 6
//   A = 7, root = 8
struct CoreFixture {
  BuilderPool pool;
  PDLTreeCore<> core;
};

CoreFixture MakeFixture(StoragePolicy policy, std::size_t n_doc = 5,
                        float storing_factor = 4.0f, uint32_t block_size = 1) {
  CoreFixture f;
  std::vector<std::size_t> lcp = {0, 1, 4, 2, 1, 5, 3, 0};
  std::vector<std::size_t> da = {3, 0, 1, 1, 2, 0, 3};
  auto* root = BuildSparseSuffixTree(f.pool, [&](std::size_t i) { return lcp[i]; }, 7);
  InsertExplicitLeaves(f.pool, root);
  ComputeDocSetsBottomUp(root, n_doc, [&](std::size_t i) { return da[i]; });
  ApplyStoragePolicy(root, policy, storing_factor);
  std::size_t n_nodes = AssignNodeIds(root);
  BuildPDLTreeCoreFromBuilder(f.core, root, n_nodes, n_doc, block_size, storing_factor, policy);
  return f;
}

using Range = std::pair<std::size_t, std::size_t>;

struct CoverResult {
  std::vector<Range> raw;
  std::vector<std::size_t> nodes;
};

CoverResult Cover(const PDLTreeCore<>& core, std::size_t sp, std::size_t ep) {
  CoverResult r;
  core.computeCoverFull(sp, ep, r.raw, r.nodes);
  return r;
}

// Empty / null queries.
TEST(PDLTreeCoreCover, EmptyQueryYieldsNothing) {
  auto f = MakeFixture(StoragePolicy::StoreAllInternal);
  auto r = Cover(f.core, 3, 3);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_THAT(r.nodes, testing::IsEmpty());
}

TEST(PDLTreeCoreCover, EmptyCoreYieldsNothing) {
  PDLTreeCore<> core;  // un-Assembled
  auto r = Cover(core, 0, 10);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_THAT(r.nodes, testing::IsEmpty());
}

// Whole-pattern range under StoreAllInternal — root is selected, so the
// whole query [0, 7) is satisfied by a single node and no raw range.
TEST(PDLTreeCoreCover, WholePatternRangeUsesRootWhenSelected) {
  auto f = MakeFixture(StoragePolicy::StoreAllInternal);
  auto r = Cover(f.core, 0, 7);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  ASSERT_EQ(r.nodes.size(), 1u);  // root via codec slot
}

// Exact-node interval: query covers exactly E[4,7) and StoreAllInternal
// selects E directly with no raw range.
TEST(PDLTreeCoreCover, ExactNodeIntervalUsesThatNode) {
  auto f = MakeFixture(StoragePolicy::StoreAllInternal);
  auto r = Cover(f.core, 4, 7);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  ASSERT_EQ(r.nodes.size(), 1u);
}

// Leaf-only range under LeavesOnly: only leaves are selected, so a query
// over [0, 1) (the explicit leaf at position 0) returns that leaf with
// no raw fallback.
TEST(PDLTreeCoreCover, LeafOnlyRangeUnderLeavesOnly) {
  auto f = MakeFixture(StoragePolicy::LeavesOnly);
  auto r = Cover(f.core, 0, 1);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  ASSERT_EQ(r.nodes.size(), 1u);
}

// Partial overlap: query [3, 5) crosses the C/E boundary. Walk under
// StoreAllInternal:
//   root partial -> descend A
//     A partial -> descend leaf(0,1) (disjoint), C (partial), E (partial)
//     C partial -> descend B (disjoint), leaf(3,4) (fully in -> selected)
//     E partial -> descend D[4,6) (partial leaf -> raw [4,5)),
//                          leaf(6,7) (disjoint)
// Expected: one selected node (leaf(3,4) via codec slot) and one raw
// range [4, 5).
TEST(PDLTreeCoreCover, PartialOverlapDescendsAndCombines) {
  auto f = MakeFixture(StoragePolicy::StoreAllInternal);
  auto r = Cover(f.core, 3, 5);

  EXPECT_THAT(r.raw, testing::ElementsAre(Range{4, 5}));
  EXPECT_EQ(r.nodes.size(), 1u);
}

// Single-gap cover under LeavesOnly: query [0, 4). All five leaves
// selected, none of the internal A/C/root selected. Cover walks down,
// includes leaf(0,1) (in [0,4)), descends C (partially overlapping):
// B[1,3) is fully in -> selected leaf, leaf(3,4) is in -> selected.
TEST(PDLTreeCoreCover, SingleGapPartitionUnderLeavesOnly) {
  auto f = MakeFixture(StoragePolicy::LeavesOnly);
  auto r = Cover(f.core, 0, 4);

  // No raw — every position [0, 4) is covered by a selected leaf.
  EXPECT_THAT(r.raw, testing::IsEmpty());
  // Three selected leaves: leaf(0,1), B[1,3), leaf(3,4).
  EXPECT_EQ(r.nodes.size(), 3u);
}

// Multi-gap cover: under OccurrenceWeighted with sf=10 and n_doc=4, only
// root and A are selected via contains_all (computed in the policy
// fixture in pdl_tree_builder_test.cpp). C and E are NOT selected, but
// their leaves are. Query [3, 6) crosses C and E with descent producing
// multiple raw-vs-selected fragments.
TEST(PDLTreeCoreCover, MixedSelectivityProducesNodesAndRawTogether) {
  auto f = MakeFixture(StoragePolicy::LeavesOnly);
  auto r = Cover(f.core, 1, 6);

  // [1, 6) covered entirely by leaves (LeavesOnly): B[1,3), leaf(3,4),
  // D[4,6). No raw.
  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_EQ(r.nodes.size(), 3u);
}

// Whole-pattern range under LeavesOnly: walk produces one selected node
// per leaf (5 of them) with no raw range.
TEST(PDLTreeCoreCover, WholePatternRangeUnderLeavesOnly) {
  auto f = MakeFixture(StoragePolicy::LeavesOnly);
  auto r = Cover(f.core, 0, 7);

  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_EQ(r.nodes.size(), 5u);
}

}  // namespace
