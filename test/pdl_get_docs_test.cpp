//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/11/26.
//
// PDLRawRangePolicy interface tests (Task 19). The wrapper is a thin
// adapter around the rmq::GetDoc* family; here we use a hand-written
// stand-in so the contract can be exercised without spinning up a full
// PDL index. Tasks 20-22 add tests against the real DA / GCDA / DGCDA
// backings.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "dret/pdl/get_docs.h"

namespace {

using dret::pdl::PDLRawRangePolicy;

// Mock raw-range source matching the rmq::GetDoc* shape: operator()(i)
// returns one doc, operator()(b, e, report) iterates a half-open range.
struct MockGetDoc {
  std::vector<std::size_t> da;

  std::size_t operator()(std::size_t i) const { return da.at(i); }

  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& r) const {
    for (std::size_t i = b; i < e; ++i) r(da.at(i));
  }
};

PDLRawRangePolicy<MockGetDoc> MakePolicy(std::vector<std::size_t> da) {
  return PDLRawRangePolicy<MockGetDoc>{MockGetDoc{std::move(da)}};
}

TEST(PDLRawRangePolicy, GetDocSinglePosition) {
  auto policy = MakePolicy({3, 0, 1, 1, 2});
  EXPECT_EQ(policy.getDoc(0), 3u);
  EXPECT_EQ(policy.getDoc(2), 1u);
  EXPECT_EQ(policy.getDoc(4), 2u);
}

TEST(PDLRawRangePolicy, GetDocsHalfOpenRange) {
  auto policy = MakePolicy({3, 0, 1, 1, 2});
  std::vector<std::size_t> seen;
  auto report = [&seen](std::size_t d) { seen.push_back(d); };
  policy.getDocs(/*sp=*/1, /*ep=*/4, report);
  EXPECT_THAT(seen, testing::ElementsAre(0u, 1u, 1u));
}

TEST(PDLRawRangePolicy, GetDocsEmptyRangeReportsNothing) {
  auto policy = MakePolicy({3, 0, 1, 1, 2});
  std::vector<std::size_t> seen;
  auto report = [&seen](std::size_t d) { seen.push_back(d); };
  policy.getDocs(/*sp=*/2, /*ep=*/2, report);
  EXPECT_THAT(seen, testing::IsEmpty());
}

TEST(PDLRawRangePolicy, GetDocsFullRange) {
  auto policy = MakePolicy({3, 0, 1, 1, 2});
  std::vector<std::size_t> seen;
  auto report = [&seen](std::size_t d) { seen.push_back(d); };
  policy.getDocs(/*sp=*/0, /*ep=*/5, report);
  EXPECT_THAT(seen, testing::ElementsAre(3u, 0u, 1u, 1u, 2u));
}

TEST(PDLRawRangePolicy, InnerAccessExposesWrappedObject) {
  auto policy = MakePolicy({7, 7, 7});
  EXPECT_EQ(policy.inner().da.size(), 3u);
  // Mutate via inner() to confirm the adapter holds, not copies.
  policy.inner().da[0] = 9;
  EXPECT_EQ(policy.getDoc(0), 9u);
}

TEST(PDLRawRangePolicy, DefaultConstructibleHoldsDefaultInner) {
  PDLRawRangePolicy<MockGetDoc> policy;
  EXPECT_TRUE(policy.inner().da.empty());
}

// Compile-checks for Tasks 20 / 21 / 22: the three rmq::GetDoc* aliases
// must instantiate at default template arguments. Real cross-type
// equivalence ("DA, GCDA, DGCDA report identical doc ids over identical
// ranges") is the subject of Task 38; here we only pin that the
// adapter wiring compiles.
TEST(PDLRawRangePolicy, AliasesCompileAtDefaults) {
  using DA = dret::pdl::PDLGetDocsDA<>;
  using SLP = dret::pdl::PDLGetDocsSLP<>;
  using DSLP = dret::pdl::PDLGetDocsDSLP<>;

  static_assert(sizeof(DA) > 0);
  static_assert(sizeof(SLP) > 0);
  static_assert(sizeof(DSLP) > 0);

  // Default-construct each alias; at this point inner() returns a
  // default-constructed rmq::GetDoc* (no cache loaded), so we can't
  // safely call getDocs(...) -- that's what Task 38's integration tests
  // will exercise. Just confirm construction doesn't fault.
  DA da;
  SLP slp;
  DSLP dslp;
  EXPECT_EQ(&da.inner(), &da.inner());
  EXPECT_EQ(&slp.inner(), &slp.inner());
  EXPECT_EQ(&dslp.inner(), &dslp.inner());
}

}  // namespace
