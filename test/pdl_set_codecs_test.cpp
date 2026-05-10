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

}  // namespace
