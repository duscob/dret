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

}  // namespace
