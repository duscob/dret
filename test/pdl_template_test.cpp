//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/11/26.
//
// Task 18: template-instantiation tests for the bitvector axis. Each
// typed test parameterises a class (PDLTreeCore or BCCodec) on
// TBitvector ∈ {sdsl::sd_vector<>, sdsl::bit_vector} and asserts that
// the configured instantiation compiles AND produces identical query
// results to the default. Catches both compile-time regressions in
// templates and runtime regressions in the rank/select rebind path.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <sstream>
#include <utility>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/sd_vector.hpp>

#include "dret/pdl/build_pdl_core.h"
#include "dret/pdl/set_codecs.h"
#include "dret/pdl/storage_policy.h"
#include "dret/pdl/tree_builder.h"
#include "dret/pdl/tree_core.h"

namespace {

using dret::pdl::ApplyStoragePolicy;
using dret::pdl::AssignNodeIds;
using dret::pdl::BCCodec;
using dret::pdl::BuildPDLTreeCoreFromBuilder;
using dret::pdl::BuildSparseSuffixTree;
using dret::pdl::BuilderPool;
using dret::pdl::ComputeDocSetsBottomUp;
using dret::pdl::InsertExplicitLeaves;
using dret::pdl::PDLTreeCore;
using dret::pdl::PlainCodec;
using dret::pdl::StoragePolicy;
using dret::pdl::StoredSet;

//===========================================================================
// PDLTreeCore<TBitvector> typed tests
//===========================================================================

template <typename TBitvector>
struct PDLTreeCoreBitvectorTraits {
  using Bitvector = TBitvector;
  using Core = PDLTreeCore<TBitvector>;
};

template <typename TTraits>
class PDLTreeCoreBitvectorTypedTests : public ::testing::Test {};

using PDLTreeCoreBitvectorTraitsList = ::testing::Types<
    PDLTreeCoreBitvectorTraits<sdsl::sd_vector<>>,
    PDLTreeCoreBitvectorTraits<sdsl::bit_vector>>;

TYPED_TEST_SUITE(PDLTreeCoreBitvectorTypedTests, PDLTreeCoreBitvectorTraitsList);

// Build the canonical PDL pipeline fixture (lcp = [0,1,4,2,1,5,3,0],
// n=7, DA = {3,0,1,1,2,0,3}) and bridge into a PDLTreeCore<TBitvector>.
template <typename TCore>
TCore BuildFixture(BuilderPool& pool, StoragePolicy policy) {
  std::vector<std::size_t> lcp = {0, 1, 4, 2, 1, 5, 3, 0};
  std::vector<std::size_t> da = {3, 0, 1, 1, 2, 0, 3};
  auto* root = BuildSparseSuffixTree(pool, [&](std::size_t i) { return lcp[i]; }, 7);
  InsertExplicitLeaves(pool, root);
  ComputeDocSetsBottomUp(root, /*n_doc=*/5, [&](std::size_t i) { return da[i]; });
  ApplyStoragePolicy(root, policy, 4.0f);
  std::size_t n_nodes = AssignNodeIds(root);

  TCore core;
  BuildPDLTreeCoreFromBuilder(core, root, n_nodes, /*n_doc=*/5,
                              /*block_size=*/1, /*storing_factor=*/4.0f, policy);
  return core;
}

using Range = std::pair<std::size_t, std::size_t>;
struct CoverProbe {
  std::vector<Range> raw;
  std::vector<std::size_t> nodes;
  bool operator==(const CoverProbe&) const = default;
};

template <typename TCore>
CoverProbe Cover(const TCore& core, std::size_t sp, std::size_t ep) {
  CoverProbe p;
  core.computeCoverFull(sp, ep, p.raw, p.nodes);
  return p;
}

template <typename TCore>
std::vector<CoverProbe> ProbeAllQueries(const TCore& core) {
  std::vector<CoverProbe> out;
  for (std::size_t sp = 0; sp <= 7; ++sp) {
    for (std::size_t ep = sp; ep <= 7; ++ep) {
      out.push_back(Cover(core, sp, ep));
    }
  }
  return out;
}

// Reference (sd_vector default) probe: build once and reuse for parity
// checks across bitvector variants. Exposed as a function to keep the
// test bodies self-contained.
std::vector<CoverProbe> ReferenceProbeStoreAllInternal() {
  BuilderPool pool;
  auto core = BuildFixture<PDLTreeCore<>>(pool, StoragePolicy::StoreAllInternal);
  return ProbeAllQueries(core);
}

std::vector<CoverProbe> ReferenceProbeLeavesOnly() {
  BuilderPool pool;
  auto core = BuildFixture<PDLTreeCore<>>(pool, StoragePolicy::LeavesOnly);
  return ProbeAllQueries(core);
}

TYPED_TEST(PDLTreeCoreBitvectorTypedTests, CompilesAndProducesReferenceCoverStoreAllInternal) {
  using Core = typename TypeParam::Core;
  BuilderPool pool;
  auto core = BuildFixture<Core>(pool, StoragePolicy::StoreAllInternal);

  EXPECT_EQ(ProbeAllQueries(core), ReferenceProbeStoreAllInternal());
}

TYPED_TEST(PDLTreeCoreBitvectorTypedTests, CompilesAndProducesReferenceCoverLeavesOnly) {
  using Core = typename TypeParam::Core;
  BuilderPool pool;
  auto core = BuildFixture<Core>(pool, StoragePolicy::LeavesOnly);

  EXPECT_EQ(ProbeAllQueries(core), ReferenceProbeLeavesOnly());
}

TYPED_TEST(PDLTreeCoreBitvectorTypedTests, SerializeLoadRoundTrip) {
  using Core = typename TypeParam::Core;
  BuilderPool pool;
  auto core = BuildFixture<Core>(pool, StoragePolicy::StoreAllInternal);
  auto before = ProbeAllQueries(core);

  std::stringstream ss;
  std::size_t bytes = core.serialize(ss);
  EXPECT_GT(bytes, 0u);

  Core reloaded;
  reloaded.load(ss);

  EXPECT_EQ(ProbeAllQueries(reloaded), before);
}

//===========================================================================
// BCCodec<TBitvector> typed tests
//===========================================================================

template <typename TBitvector>
struct BCCodecBitvectorTraits {
  using Bitvector = TBitvector;
  using Codec = BCCodec<TBitvector>;
};

template <typename TTraits>
class BCCodecBitvectorTypedTests : public ::testing::Test {};

using BCCodecBitvectorTraitsList = ::testing::Types<
    BCCodecBitvectorTraits<sdsl::sd_vector<>>,
    BCCodecBitvectorTraits<sdsl::bit_vector>>;

TYPED_TEST_SUITE(BCCodecBitvectorTypedTests, BCCodecBitvectorTraitsList);

template <typename TCodec>
std::vector<std::size_t> ExpandSorted(const TCodec& codec, std::size_t slot,
                                      std::size_t n_doc) {
  std::vector<std::size_t> out;
  codec.Expand(slot, n_doc, [&out](std::size_t d) { out.push_back(d); });
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Same fixture shape BCCodec parity tests use; a 9-slot mix exercising
// every shape PDLBCCodec.ExpandRepeatedPatternsExerciseRules covers
// (but we don't care here whether vnmextract finds bicliques — only
// that the typed instantiation matches the default).
const std::vector<StoredSet>& ParityFixture() {
  static const std::vector<StoredSet> sets = {
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1}},
      {true,  {}},
      {false, {2, 3}},
      {false, {1, 2}},
      {false, {0, 3}},
  };
  return sets;
}

TYPED_TEST(BCCodecBitvectorTypedTests, CompilesAndProducesReferenceExpansion) {
  using Codec = typename TypeParam::Codec;
  const auto& sets = ParityFixture();
  const std::size_t n_doc = 4;

  // Reference codec (sd_vector default).
  BCCodec<> reference;
  reference.Build(sets.size(), [&sets](std::size_t s) { return sets[s]; }, n_doc);

  // Variant under test.
  Codec variant;
  variant.Build(sets.size(), [&sets](std::size_t s) { return sets[s]; }, n_doc);

  ASSERT_EQ(variant.n_slots(), reference.n_slots());
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_EQ(ExpandSorted(variant, s, n_doc), ExpandSorted(reference, s, n_doc));
  }
}

TYPED_TEST(BCCodecBitvectorTypedTests, SerializeLoadRoundTrip) {
  using Codec = typename TypeParam::Codec;
  const auto& sets = ParityFixture();
  const std::size_t n_doc = 4;

  Codec codec;
  codec.Build(sets.size(), [&sets](std::size_t s) { return sets[s]; }, n_doc);

  std::stringstream ss;
  std::size_t bytes = codec.serialize(ss);
  EXPECT_GT(bytes, 0u);

  Codec reloaded;
  reloaded.load(ss);

  EXPECT_EQ(reloaded.n_slots(), codec.n_slots());
  EXPECT_EQ(reloaded.n_rules(), codec.n_rules());
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_EQ(ExpandSorted(reloaded, s, n_doc), ExpandSorted(codec, s, n_doc));
  }
}

}  // namespace
