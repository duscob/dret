//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// Direct PDLTreeCore::computeCoverFull tests on hand-assembled tiny
// trees, Task 35 of docs/pdl_indexes_tasks.md. Unlike
// pdl_tree_core_test.cpp (which builds the core via the full builder
// pipeline), these tests bypass the builder and call
// PDLTreeCore::Assemble directly with hand-coded int/bit vectors so
// the expected raw ranges and selected node ids are obvious from the
// fixture itself — i.e., cover behavior is validated independently of
// any tree-builder, codec, or pattern-counting correctness.
//
// Acceptance criteria (Task 35): exact interval, partial interval,
// leaf-only, and full-range covers match expected outputs.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <utility>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

#include "dret/pdl/storage_policy.h"
#include "dret/pdl/tree_core.h"

namespace {

using dret::pdl::PDLTreeCore;
using dret::pdl::StoragePolicy;

using Range = std::pair<std::size_t, std::size_t>;

struct CoverResult {
  std::vector<Range> raw;
  std::vector<std::size_t> nodes;
};

// Hand-build a 3-node tree:
//
//   root [0, 4)  id = 2
//     leaf [0, 2)  id = 0
//     leaf [2, 4)  id = 1
//
// n_nodes_ = 3, kSentinel = 3. AssignNodeIds-equivalent post-order
// gives leaf0 = 0, leaf1 = 1, root = 2 (root is always n_nodes_ - 1).
//
// `selected_marker[i] = 1` means node i has a stored set; codec slot
// is selected_rank_(id) (the rank-1 of id in selected_marker_).
//
// selected_pattern controls which nodes are selected — three of these
// fixtures correspond to the three storage policies (OccurrenceWeighted
// example, StoreAllInternal, LeavesOnly), but here policy is just
// metadata; the selected_marker is what the cover algorithm actually
// reads.
PDLTreeCore<> MakeTree(std::vector<int> selected_pattern) {
  const std::size_t kN = 3;
  const std::size_t kSent = kN;

  auto make_iv = [](std::initializer_list<std::size_t> vals) {
    sdsl::int_vector<> iv(vals.size(), 0, 64);
    std::size_t i = 0;
    for (auto v : vals) iv[i++] = v;
    return iv;
  };

  auto node_starts  = make_iv({0, 2, 0});
  auto node_ends    = make_iv({2, 4, 4});
  auto first_child  = make_iv({kSent, kSent, 0});  // only root has children
  auto next_sibling = make_iv({1, kSent, kSent});  // leaf0 -> leaf1

  sdsl::bit_vector marker_bv(kN, 0);
  for (std::size_t i = 0; i < kN && i < selected_pattern.size(); ++i) {
    marker_bv[i] = selected_pattern[i] ? 1 : 0;
  }
  sdsl::sd_vector<> marker(marker_bv);

  PDLTreeCore<> core;
  core.Assemble(std::move(node_starts), std::move(node_ends),
                std::move(first_child), std::move(next_sibling),
                std::move(marker), {},  // NullCodec — cover doesn't read it.
                /*n_doc=*/2,
                /*block_size=*/1, /*storing_factor=*/1.0f,
                StoragePolicy::OccurrenceWeighted);
  return core;
}

CoverResult Cover(const PDLTreeCore<>& core, std::size_t sp, std::size_t ep) {
  CoverResult r;
  core.computeCoverFull(sp, ep, r.raw, r.nodes);
  return r;
}

// Full-range cover — all three nodes selected so the root absorbs the
// whole query in one node hit.
TEST(PDLTreeCoverDirect, FullRangeRootSelected) {
  auto core = MakeTree({1, 1, 1});
  auto r = Cover(core, 0, 4);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  // Codec slot for root = selected_rank_(2) = count of 1-bits in [0, 2) = 2.
  EXPECT_THAT(r.nodes, testing::ElementsAre(2u));
}

// Exact-interval cover — query exactly matches a single selected node.
TEST(PDLTreeCoverDirect, ExactIntervalUsesThatNode) {
  auto core = MakeTree({1, 1, 1});
  auto r = Cover(core, 0, 2);  // exactly leaf0
  EXPECT_THAT(r.raw, testing::IsEmpty());
  // Codec slot for leaf0 = selected_rank_(0) = 0.
  EXPECT_THAT(r.nodes, testing::ElementsAre(0u));
}

// Leaf-only cover — root NOT selected, both leaves selected. The full
// query [0, 4) descends from root and finds two selected leaves.
TEST(PDLTreeCoverDirect, LeafOnlyFullRangeDescends) {
  auto core = MakeTree({1, 1, 0});
  auto r = Cover(core, 0, 4);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  // Codec slots: leaf0 -> rank_(0) = 0; leaf1 -> rank_(1) = 1.
  EXPECT_THAT(r.nodes, testing::UnorderedElementsAre(0u, 1u));
}

// Partial-interval cover — root NOT selected, query [0, 3) crosses the
// leaf0/leaf1 boundary. leaf0 fully inside -> selected slot; leaf1
// partially overlaps -> raw [2, 3) (because it's a leaf).
TEST(PDLTreeCoverDirect, PartialIntervalMixesNodeAndRaw) {
  auto core = MakeTree({1, 1, 0});
  auto r = Cover(core, 0, 3);
  EXPECT_THAT(r.raw, testing::ElementsAre(Range{2, 3}));
  EXPECT_THAT(r.nodes, testing::ElementsAre(0u));
}

// Partial-interval where neither leaf is selected — both partial leaves
// should fall back to raw ranges, no node hits.
TEST(PDLTreeCoverDirect, PartialIntervalAllRawWhenLeavesUnselected) {
  auto core = MakeTree({0, 0, 0});
  auto r = Cover(core, 1, 3);
  // Walk: root partial -> descend; leaf0 [0,2) partial -> raw [1,2);
  // leaf1 [2,4) partial -> raw [2,3). Order is stack-LIFO (leaf1
  // popped before leaf0), so assert order-agnostically.
  EXPECT_THAT(r.raw, testing::UnorderedElementsAre(Range{1, 2}, Range{2, 3}));
  EXPECT_THAT(r.nodes, testing::IsEmpty());
}

// Empty query — sp == ep returns nothing regardless of marker pattern.
TEST(PDLTreeCoverDirect, EmptyQueryYieldsNothing) {
  auto core = MakeTree({1, 1, 1});
  auto r = Cover(core, 2, 2);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_THAT(r.nodes, testing::IsEmpty());
}

// Disjoint query — entirely outside the tree extent.
TEST(PDLTreeCoverDirect, DisjointQueryYieldsNothing) {
  auto core = MakeTree({1, 1, 1});
  auto r = Cover(core, 10, 20);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_THAT(r.nodes, testing::IsEmpty());
}

// Empty (un-Assembled) core — cover must short-circuit instead of
// indexing into empty arrays.
TEST(PDLTreeCoverDirect, EmptyCoreYieldsNothing) {
  PDLTreeCore<> core;
  auto r = Cover(core, 0, 4);
  EXPECT_THAT(r.raw, testing::IsEmpty());
  EXPECT_THAT(r.nodes, testing::IsEmpty());
}

}  // namespace
