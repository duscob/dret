//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// PDL stored-set codec tests. Task 14 lands the SetCodec concept and a
// DummyCodec reference impl; Tasks 15-17 extend this file with the real
// Plain / RP / BC codecs. The include below forces the in-header
// static_assert(SetCodec<DummyCodec>) to fire — that's the compile-only
// acceptance check from the task description.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include "dret/pdl/set_codecs.h"

namespace {

using dret::pdl::BCCodec;
using dret::pdl::DummyCodec;
using dret::pdl::ExpandAllDoc;
using dret::pdl::PlainCodec;
using dret::pdl::RPCodec;
using dret::pdl::SetCodec;
using dret::pdl::StoredSet;

// Belt-and-braces: re-assert the concept holds for DummyCodec inside
// translation-unit context. If the SetCodec definition ever drifts,
// either this assert or the in-header one will catch it.
static_assert(SetCodec<DummyCodec>);

TEST(PDLDummyCodec, ExpandEmitsAllDocsZeroBased) {
  DummyCodec codec;
  std::vector<std::size_t> seen;
  codec.Expand(/*slot=*/0, /*n_doc=*/4,
               [&seen](std::size_t d) { seen.push_back(d); });
  EXPECT_THAT(seen, testing::ElementsAre(0u, 1u, 2u, 3u));
}

TEST(PDLDummyCodec, ExpandWithZeroDocsEmitsNothing) {
  DummyCodec codec;
  std::vector<std::size_t> seen;
  codec.Expand(/*slot=*/0, /*n_doc=*/0,
               [&seen](std::size_t d) { seen.push_back(d); });
  EXPECT_THAT(seen, testing::IsEmpty());
}

TEST(PDLDummyCodec, BuildAcceptsCallable) {
  DummyCodec codec;
  // Build is a no-op on DummyCodec; just verify the call shape compiles
  // for a slot-indexed StoredSet source — the signature future codecs
  // will implement.
  codec.Build(/*n_slots=*/3, [](std::size_t slot) {
    StoredSet s;
    s.contains_all = (slot == 0);
    if (!s.contains_all) s.docs = {slot, slot + 1};
    return s;
  });
}

TEST(PDLDummyCodec, SerializeLoadRoundTrip) {
  DummyCodec codec;
  std::stringstream ss;
  std::size_t bytes = codec.serialize(ss);
  EXPECT_EQ(bytes, 0u);  // dummy stores nothing

  DummyCodec reloaded;
  reloaded.load(ss);

  std::vector<std::size_t> seen;
  reloaded.Expand(0, 2, [&seen](std::size_t d) { seen.push_back(d); });
  EXPECT_THAT(seen, testing::ElementsAre(0u, 1u));
}

TEST(PDLExpandAllDoc, EmitsRange) {
  std::vector<std::size_t> seen;
  ExpandAllDoc(5, [&seen](std::size_t d) { seen.push_back(d); });
  EXPECT_THAT(seen, testing::ElementsAre(0u, 1u, 2u, 3u, 4u));
}

// PlainCodec tests (Task 15). Acceptance criterion: codec expansion
// returns exact input sets and expands all-doc sentinel to 0..n_doc-1.

static_assert(SetCodec<PlainCodec<>>);

// Build helper: converts an in-memory list of StoredSet records into a
// slot-indexed source the codec's Build can consume.
auto SetSource(const std::vector<StoredSet>& sets) {
  return [&sets](std::size_t slot) { return sets[slot]; };
}

std::vector<std::size_t> ExpandAt(const PlainCodec<>& codec, std::size_t slot,
                                  std::size_t n_doc) {
  std::vector<std::size_t> out;
  codec.Expand(slot, n_doc, [&out](std::size_t d) { out.push_back(d); });
  return out;
}

TEST(PDLPlainCodec, ExpandSingletonSet) {
  std::vector<StoredSet> sets = {{false, {7}}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/10);

  EXPECT_EQ(codec.n_slots(), 1u);
  EXPECT_THAT(ExpandAt(codec, 0, 10), testing::ElementsAre(7u));
}

TEST(PDLPlainCodec, ExpandMultiDocSet) {
  std::vector<StoredSet> sets = {{false, {0, 2, 4, 9}}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/10);

  EXPECT_THAT(ExpandAt(codec, 0, 10), testing::ElementsAre(0u, 2u, 4u, 9u));
}

TEST(PDLPlainCodec, ExpandEmptySetEmitsNothing) {
  std::vector<StoredSet> sets = {{false, {}}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/10);

  EXPECT_THAT(ExpandAt(codec, 0, 10), testing::IsEmpty());
}

TEST(PDLPlainCodec, ExpandAllDocSentinel) {
  std::vector<StoredSet> sets = {{true, {}}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/4);

  EXPECT_THAT(ExpandAt(codec, 0, 4), testing::ElementsAre(0u, 1u, 2u, 3u));
}

// Single-element non-sentinel set whose value happens to equal n_doc-1
// (the largest legal doc id) should NOT be confused with the sentinel
// (which uses value n_doc, not n_doc-1).
TEST(PDLPlainCodec, ExpandLargestLegalSingleDocIsNotSentinel) {
  std::vector<StoredSet> sets = {{false, {9}}};  // 9 is legal when n_doc=10
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/10);

  EXPECT_THAT(ExpandAt(codec, 0, 10), testing::ElementsAre(9u));
}

// Mixed-slot fixture exercising every storage shape at once.
TEST(PDLPlainCodec, ExpandMixedFixtureAllSlots) {
  std::vector<StoredSet> sets = {
      {false, {0}},               // slot 0: singleton
      {false, {1, 2, 3}},         // slot 1: multi
      {false, {}},                // slot 2: empty
      {true,  {}},                // slot 3: all-doc sentinel
      {false, {3}},               // slot 4: largest legal singleton
  };
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/4);

  EXPECT_EQ(codec.n_slots(), 5u);
  EXPECT_THAT(ExpandAt(codec, 0, 4), testing::ElementsAre(0u));
  EXPECT_THAT(ExpandAt(codec, 1, 4), testing::ElementsAre(1u, 2u, 3u));
  EXPECT_THAT(ExpandAt(codec, 2, 4), testing::IsEmpty());
  EXPECT_THAT(ExpandAt(codec, 3, 4), testing::ElementsAre(0u, 1u, 2u, 3u));
  EXPECT_THAT(ExpandAt(codec, 4, 4), testing::ElementsAre(3u));
}

TEST(PDLPlainCodec, SerializeLoadRoundTripPreservesExpansions) {
  std::vector<StoredSet> sets = {
      {false, {0, 4}},
      {true,  {}},
      {false, {2}},
  };
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/5);

  std::stringstream ss;
  std::size_t bytes = codec.serialize(ss);
  EXPECT_GT(bytes, 0u);

  PlainCodec<> reloaded;
  reloaded.load(ss);

  EXPECT_EQ(reloaded.n_slots(), codec.n_slots());
  EXPECT_THAT(ExpandAt(reloaded, 0, 5), testing::ElementsAre(0u, 4u));
  EXPECT_THAT(ExpandAt(reloaded, 1, 5), testing::ElementsAre(0u, 1u, 2u, 3u, 4u));
  EXPECT_THAT(ExpandAt(reloaded, 2, 5), testing::ElementsAre(2u));
}

TEST(PDLPlainCodec, GetSizeReportNonzero) {
  std::vector<StoredSet> sets = {{false, {1, 2, 3}}, {true, {}}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/5);

  auto report = codec.GetSizeReport();
  ASSERT_FALSE(report.empty());
  for (const auto& f : report) {
    EXPECT_GT(f.bytes, 0u) << "field " << f.name;
  }
}

// RPCodec tests (Task 16). Acceptance: expansion matches PlainCodec for
// the same selected sets, including nested rules and the all-doc
// sentinel.

static_assert(SetCodec<RPCodec<>>);

std::vector<std::size_t> RPExpandAt(const RPCodec<>& codec, std::size_t slot,
                                    std::size_t n_doc) {
  std::vector<std::size_t> out;
  codec.Expand(slot, n_doc, [&out](std::size_t d) { out.push_back(d); });
  return out;
}

// Convenience: build both codecs from the same set source and confirm
// each slot's expansion matches.
void ExpectPlainRPParity(const std::vector<StoredSet>& sets, std::size_t n_doc) {
  PlainCodec<> plain;
  RPCodec<> rp;
  plain.Build(sets.size(), SetSource(sets), n_doc);
  rp.Build(sets.size(), SetSource(sets), n_doc);

  ASSERT_EQ(plain.n_slots(), rp.n_slots());
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_EQ(RPExpandAt(rp, s, n_doc), ExpandAt(plain, s, n_doc));
  }
}

TEST(PDLRPCodec, ExpandSingletonSet) {
  ExpectPlainRPParity({{false, {7}}}, /*n_doc=*/10);
}

TEST(PDLRPCodec, ExpandMultiDocSet) {
  ExpectPlainRPParity({{false, {0, 2, 4, 9}}}, /*n_doc=*/10);
}

TEST(PDLRPCodec, ExpandAllDocSentinel) {
  ExpectPlainRPParity({{false, {0, 1}}, {true, {}}, {false, {2}}}, /*n_doc=*/4);
}

TEST(PDLRPCodec, ExpandLargestLegalSingleDocIsNotSentinel) {
  // The largest legal doc id is n_doc-1 = 9; the sentinel is 10.
  ExpectPlainRPParity({{false, {9}}}, /*n_doc=*/10);
}

// Nested-rule fixture: many slots with overlapping subsets so RePair
// has material to discover repeated pairs and build multi-level rules.
TEST(PDLRPCodec, ExpandRepeatedPatternsExerciseNestedRules) {
  std::vector<StoredSet> sets = {
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1}},
      {false, {2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {1, 2}},
      {false, {0, 3}},
  };
  ExpectPlainRPParity(sets, /*n_doc=*/4);
}

TEST(PDLRPCodec, ExpandMixedFixtureMatchesPlain) {
  // Same shape as PlainCodec's mixed fixture; tests parity end-to-end.
  std::vector<StoredSet> sets = {
      {false, {0}},
      {false, {1, 2, 3}},
      {true,  {}},
      {false, {3}},
      {false, {0, 1, 2, 3}},
  };
  ExpectPlainRPParity(sets, /*n_doc=*/4);
}

TEST(PDLRPCodec, SerializeLoadRoundTripPreservesExpansions) {
  std::vector<StoredSet> sets = {
      {false, {0, 4}},
      {true,  {}},
      {false, {2}},
      {false, {0, 1, 2, 3, 4}},
  };
  RPCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/5);

  std::stringstream ss;
  std::size_t bytes = codec.serialize(ss);
  EXPECT_GT(bytes, 0u);

  RPCodec<> reloaded;
  reloaded.load(ss);

  EXPECT_EQ(reloaded.n_slots(), codec.n_slots());
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_EQ(RPExpandAt(reloaded, s, 5), RPExpandAt(codec, s, 5));
  }
}

TEST(PDLRPCodec, GetSizeReportNonzero) {
  std::vector<StoredSet> sets = {{false, {1, 2, 3}}, {true, {}}};
  RPCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/5);

  auto report = codec.GetSizeReport();
  ASSERT_FALSE(report.empty());
  for (const auto& f : report) {
    EXPECT_GT(f.bytes, 0u) << "field " << f.name;
  }
}

// BCCodec tests (Tasks 17.B, 17.C, 17a). Acceptance: expansion matches
// PlainCodec for the same selected sets, including the all-doc sentinel
// and inputs where vnmextract finds bicliques to compress.

static_assert(SetCodec<BCCodec<>>);

std::vector<std::size_t> BCExpandAt(const BCCodec<>& codec, std::size_t slot,
                                    std::size_t n_doc) {
  std::vector<std::size_t> out;
  codec.Expand(slot, n_doc, [&out](std::size_t d) { out.push_back(d); });
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Compare per-slot Expand of PlainCodec and BCCodec (both deduped + sorted).
void ExpectPlainBCParity(const std::vector<StoredSet>& sets, std::size_t n_doc) {
  PlainCodec<> plain;
  BCCodec<> bc;
  plain.Build(sets.size(), SetSource(sets), n_doc);
  bc.Build(sets.size(), SetSource(sets), n_doc);

  ASSERT_EQ(plain.n_slots(), bc.n_slots());
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    auto plain_out = ExpandAt(plain, s, n_doc);
    std::sort(plain_out.begin(), plain_out.end());
    plain_out.erase(std::unique(plain_out.begin(), plain_out.end()), plain_out.end());
    EXPECT_EQ(BCExpandAt(bc, s, n_doc), plain_out);
  }
}

TEST(PDLBCCodec, ExpandSingletonSet) {
  ExpectPlainBCParity({{false, {7}}}, /*n_doc=*/10);
}

TEST(PDLBCCodec, ExpandMultiDocSet) {
  ExpectPlainBCParity({{false, {0, 2, 4, 9}}}, /*n_doc=*/10);
}

TEST(PDLBCCodec, ExpandEmptySetEmitsNothing) {
  std::vector<StoredSet> sets = {{false, {}}};
  BCCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/10);
  EXPECT_THAT(BCExpandAt(codec, 0, 10), testing::IsEmpty());
}

TEST(PDLBCCodec, ExpandAllDocSentinel) {
  ExpectPlainBCParity({{false, {0, 1}}, {true, {}}, {false, {2}}}, /*n_doc=*/4);
}

TEST(PDLBCCodec, ExpandLargestLegalSingleDocIsNotSentinel) {
  ExpectPlainBCParity({{false, {9}}}, /*n_doc=*/10);
}

// Repeated-pattern fixture so vnmextract has dense subsets to find. The
// parity check works whether or not bicliques are extracted — Expand
// must always agree with PlainCodec.
TEST(PDLBCCodec, ExpandRepeatedPatternsExerciseRules) {
  std::vector<StoredSet> sets = {
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1}},
      {false, {2, 3}},
      {false, {1, 2}},
      {false, {0, 3}},
  };
  ExpectPlainBCParity(sets, /*n_doc=*/4);
}

TEST(PDLBCCodec, ExpandMixedFixtureMatchesPlain) {
  std::vector<StoredSet> sets = {
      {false, {0}},
      {false, {1, 2, 3}},
      {true,  {}},
      {false, {3}},
      {false, {0, 1, 2, 3}},
  };
  ExpectPlainBCParity(sets, /*n_doc=*/4);
}

TEST(PDLBCCodec, SerializeLoadRoundTripPreservesExpansions) {
  std::vector<StoredSet> sets = {
      {false, {0, 4}},
      {true,  {}},
      {false, {2}},
      {false, {0, 1, 2, 3, 4}},
  };
  BCCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/5);

  std::stringstream ss;
  std::size_t bytes = codec.serialize(ss);
  EXPECT_GT(bytes, 0u);

  BCCodec<> reloaded;
  reloaded.load(ss);

  EXPECT_EQ(reloaded.n_slots(), codec.n_slots());
  EXPECT_EQ(reloaded.n_rules(), codec.n_rules());
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_EQ(BCExpandAt(reloaded, s, 5), BCExpandAt(codec, s, 5));
  }
}

TEST(PDLBCCodec, GetSizeReportNonzero) {
  std::vector<StoredSet> sets = {{false, {1, 2, 3}}, {true, {}}};
  BCCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/5);

  auto report = codec.GetSizeReport();
  ASSERT_FALSE(report.empty());
}

// Task 36 additions — fill the empty/repeated/large gaps and pin the
// "expand never reports sentinel n_doc" invariant for every codec.

// --- expansion invariants helper ---
//
// Calls Expand(slot) directly (no sort/dedup) so we observe the raw
// stream the codec produces, then asserts every emitted value is
// strictly less than n_doc. This is the load-bearing acceptance
// criterion of Task 36: the all-doc sentinel must be expanded into
// 0..n_doc-1, never leaked verbatim.
template <typename TCodec>
void ExpectNoSentinelEverEmitted(const TCodec& codec, std::size_t n_slots,
                                 std::size_t n_doc) {
  for (std::size_t s = 0; s < n_slots; ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    codec.Expand(s, n_doc, [n_doc](std::size_t d) {
      EXPECT_LT(d, n_doc);
    });
  }
}

// --- repeated-set fixture for Plain (parity-style smoke for the
// pattern that already has explicit RP/BC equivalents) ---
TEST(PDLPlainCodec, ExpandRepeatedSetsAreIndependent) {
  std::vector<StoredSet> sets = {
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
      {false, {0, 1, 2, 3}},
  };
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/4);

  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_THAT(ExpandAt(codec, s, 4), testing::ElementsAre(0u, 1u, 2u, 3u));
  }
}

// --- empty-slot test for RP (gap; PlainCodec / BCCodec already have
// this). RePair upstream cannot encode a fully empty input stream, so
// the empty slot is exercised within a mixed fixture (one empty + one
// non-empty slot). The non-empty slot gives RePair material to build a
// minimal grammar; the empty slot's Expand must still emit nothing.
TEST(PDLRPCodec, ExpandEmptySlotInMixedFixtureEmitsNothing) {
  std::vector<StoredSet> sets = {
      {false, {}},          // slot 0: empty
      {false, {2, 5, 7}},   // slot 1: non-empty so RePair has input
  };
  RPCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/10);

  EXPECT_THAT(RPExpandAt(codec, 0, 10), testing::IsEmpty());
  EXPECT_THAT(RPExpandAt(codec, 1, 10), testing::ElementsAre(2u, 5u, 7u));
}

// --- large-set expansion across all three codecs ---
TEST(PDLPlainCodec, ExpandLargeSet) {
  constexpr std::size_t kNDoc = 64;
  std::vector<std::size_t> docs;
  for (std::size_t i = 0; i < kNDoc; ++i) docs.push_back(i);
  std::vector<StoredSet> sets = {{false, docs}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), kNDoc);

  auto out = ExpandAt(codec, 0, kNDoc);
  ASSERT_EQ(out.size(), kNDoc);
  for (std::size_t i = 0; i < kNDoc; ++i) EXPECT_EQ(out[i], i);
  ExpectNoSentinelEverEmitted(codec, sets.size(), kNDoc);
}

TEST(PDLRPCodec, ExpandLargeSet) {
  constexpr std::size_t kNDoc = 64;
  std::vector<std::size_t> docs;
  for (std::size_t i = 0; i < kNDoc; ++i) docs.push_back(i);
  ExpectPlainRPParity({{false, docs}}, kNDoc);

  RPCodec<> rp;
  std::vector<StoredSet> sets = {{false, docs}};
  rp.Build(sets.size(), SetSource(sets), kNDoc);
  ExpectNoSentinelEverEmitted(rp, sets.size(), kNDoc);
}

TEST(PDLBCCodec, ExpandLargeSet) {
  constexpr std::size_t kNDoc = 64;
  std::vector<std::size_t> docs;
  for (std::size_t i = 0; i < kNDoc; ++i) docs.push_back(i);
  ExpectPlainBCParity({{false, docs}}, kNDoc);

  BCCodec<> bc;
  std::vector<StoredSet> sets = {{false, docs}};
  bc.Build(sets.size(), SetSource(sets), kNDoc);
  ExpectNoSentinelEverEmitted(bc, sets.size(), kNDoc);
}

// --- "no sentinel emitted" across the all-doc-sentinel fixtures ---
TEST(PDLPlainCodec, AllDocSentinelExpandsBelowNDoc) {
  std::vector<StoredSet> sets = {{true, {}}, {false, {0, 1}}, {true, {}}};
  PlainCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/8);
  ExpectNoSentinelEverEmitted(codec, sets.size(), 8);
}

TEST(PDLRPCodec, AllDocSentinelExpandsBelowNDoc) {
  std::vector<StoredSet> sets = {{true, {}}, {false, {0, 1}}, {true, {}}};
  RPCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/8);
  ExpectNoSentinelEverEmitted(codec, sets.size(), 8);
}

TEST(PDLBCCodec, AllDocSentinelExpandsBelowNDoc) {
  std::vector<StoredSet> sets = {{true, {}}, {false, {0, 1}}, {true, {}}};
  BCCodec<> codec;
  codec.Build(sets.size(), SetSource(sets), /*n_doc=*/8);
  ExpectNoSentinelEverEmitted(codec, sets.size(), 8);
}

// --- BC-specific: dictionary (rules), codewords (blocks), borders
// must populate when vnmextract finds bicliques, and round-trip
// cleanly through serialize/load.
//
// 12 slots sharing the same 8-element doc set is dense enough for vnmextract to
// extract at least one biclique — but only under fixture-scale mining
// parameters, which this test now sets explicitly. The codec's defaults are
// drl's production values (min_bicliques=500, bcsizes=5000,500,100,50,30,15),
// which find nothing in a 12-slot fixture. That is deliberate: those values
// used to be the default for real collections too, costing >55 h on a single DA
// row of a 671 MB collection that takes 218 s at drl's settings, for 0.15% less
// space. See docs/bug_pdlbc_parameters.md.
//
// If the >0 guard ever ceases to hold (e.g. a new vnmextract version), the test
// still validates the round-trip; only that guard becomes a soft signal.
TEST(PDLBCCodec, RulesAndBlocksRoundTripWhenBicliquesExtracted) {
  constexpr std::size_t kNDoc = 16;
  std::vector<StoredSet> sets;
  // Every slot shares the high block {8..13} -- dense enough for vnmextract to
  // lift into one biclique -- plus one low doc unique to its group, which no
  // rule can cover and so stays verbatim. That partial cover is the shape that
  // matters: Expand emits the rule's docs first, then the verbatim remainder,
  // giving [8 9 10 11 12 13][0] -- correct as a set, not globally sorted.
  //
  // The earlier fixture had 12 slots sharing one identical 8-doc set, so a
  // single rule covered each slot completely, nothing was left verbatim, and
  // the expansion came out sorted by accident. It therefore could not have
  // caught the dna_010000_010 duplicate-doc defect.
  for (std::size_t i = 0; i < 12; ++i) {
    sets.push_back({false, {i % 3, 8, 9, 10, 11, 12, 13}});
  }
  BCCodec<> codec;
  codec.SetMiningParams({"1", "10,5,2", "4"});
  codec.Build(sets.size(), SetSource(sets), kNDoc);

  EXPECT_GT(codec.n_rules(), 0u)
      << "fixture intended to make vnmextract emit at least one biclique";

  // Regression: BCCodec::Expand concatenates per-rule runs with the verbatim
  // remainder, so its output is NOT globally sorted even though each run is.
  // PDLTreeCore::getDocSet feeds that to MergeSetsBinaryTreeFunctor, which
  // merges with std::set_union and requires sorted input; violating it emitted
  // 273,171 duplicate doc ids on dna_010000_010. The codec must therefore
  // declare kExpandsSorted == false so getDocSet sorts, and this pins that
  // declaration to the observed behaviour: if Expand ever starts emitting in
  // ascending order, kExpandsSorted should be revisited rather than left
  // stale.
  static_assert(!BCCodec<>::kExpandsSorted,
                "BCCodec::Expand concatenates per-rule runs, so it is not globally "
                "sorted; PDLTreeCore::getDocSet sorts on the strength of this flag");
  {
    for (std::size_t slot = 0; slot < sets.size(); ++slot) {
      std::vector<std::size_t> got;
      codec.Expand(slot, kNDoc, [&got](std::size_t d) { got.push_back(d); });
      // Deliberately not asserting sortedness: Expand concatenates runs and
      // need not be ordered. getDocSet is what must deliver sorted output, and
      // that is pinned by PDLTreeCoreGetDocSet.SortsAndUniques... using a stub
      // codec, rather than here where vnmextract's decomposition varies run to
      // run.
      // Whatever the order, the underlying set must be exactly right: the
      // dna failure had perfect sets and only duplicate emissions.
      std::vector<std::size_t> uniq(got);
      std::sort(uniq.begin(), uniq.end());
      uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
      EXPECT_EQ(uniq, sets[slot].docs) << "slot " << slot << " set mismatch";
    }
    // The set assertions above hold however vnmextract decomposed the input,
    // which is why they are the ones made here.
  }

  std::stringstream ss;
  std::size_t bytes = codec.serialize(ss);
  EXPECT_GT(bytes, 0u);

  BCCodec<> reloaded;
  reloaded.load(ss);

  EXPECT_EQ(reloaded.n_rules(), codec.n_rules());
  EXPECT_EQ(reloaded.n_slots(), codec.n_slots());

  // Per-slot expansion must round-trip; size-report must remain
  // non-empty after load (so dictionary/codeword/border vectors made
  // it through serialize).
  for (std::size_t s = 0; s < sets.size(); ++s) {
    SCOPED_TRACE(testing::Message() << "slot " << s);
    EXPECT_EQ(BCExpandAt(reloaded, s, kNDoc), BCExpandAt(codec, s, kNDoc));
  }
  ExpectNoSentinelEverEmitted(reloaded, sets.size(), kNDoc);

  auto report_after = reloaded.GetSizeReport();
  for (const auto& f : report_after) {
    if (f.name == "blocks" || f.name == "block_borders" ||
        f.name == "rules" || f.name == "rule_borders") {
      EXPECT_GT(f.bytes, 0u) << "field " << f.name;
    }
  }
}

}  // namespace
