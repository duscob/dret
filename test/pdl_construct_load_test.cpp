//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// Typed construct/load round-trip tests for every PDL index variant
// (Plain / RP / BC) crossed with every storage policy
// (OccurrenceWeighted / StoreAllInternal / LeavesOnly), Task 33 of
// docs/pdl_indexes_tasks.md. The codec is the type-list axis; the
// policy is varied at runtime inside each typed test because
// StoragePolicy is a constructor argument, not a template parameter.
//
// Acceptance: a freshly-loaded index returns the same search results
// AND the same per-field SizeReport (and total bytes) as the index it
// was built from. Search uses the DA-backed raw-range get-doc; the
// SLP/DSLP variants are covered separately by the benchmark factory
// and pdl_index_test.cpp.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/util.hpp>

#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/storage_policy.h"

#include "base_test.h"

namespace {

using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;

class DocListResultVector : public std::vector<std::size_t> {
 public:
  void operator()(std::size_t doc) { this->push_back(doc); }
  void operator()() {
    std::sort(this->begin(), this->end());
    this->erase(std::unique(this->begin(), this->end()), this->end());
  }
};

template <typename TIndex>
class DocListIdxPDLConstructLoadTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  // Same fixture used by the brute-force / smoke tests in
  // pdl_index_test.cpp; small enough that BC's vnmextract finds no
  // bicliques (verbatim path) but still exercises every codec stage.
  const std::string data_ = "TATA\1LATA\1LALA\1";

  struct PatternCase {
    std::string pattern;
    std::vector<std::size_t> docs;
  };

  const std::vector<PatternCase> patterns_ = {
      {"TAT", {0}},
      {"LAT", {1}},
      {"LAL", {2}},
      {"TA",  {0, 1}},
      {"LA",  {1, 2}},
      {"A",   {0, 1, 2}},
      {"TAL", {}},
  };

  std::vector<std::vector<std::size_t>> RunSearch(TIndex& index) {
    std::vector<std::vector<std::size_t>> all;
    all.reserve(patterns_.size());
    for (const auto& [pattern, _] : patterns_) {
      DocListResultVector result;
      index.Search(pattern, std::ref(result));
      result();
      all.emplace_back(result.begin(), result.end());
    }
    return all;
  }
};

using PDLConstructLoadTypes = ::testing::Types<
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLRP<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage>>;

TYPED_TEST_SUITE(DocListIdxPDLConstructLoadTypedTests, PDLConstructLoadTypes);

TYPED_TEST(DocListIdxPDLConstructLoadTypedTests, RoundTripPreservesSearchAndSize) {
  using Index = TypeParam;

  for (auto policy : {dret::pdl::StoragePolicy::OccurrenceWeighted,
                      dret::pdl::StoragePolicy::StoreAllInternal,
                      dret::pdl::StoragePolicy::LeavesOnly}) {
    SCOPED_TRACE(::testing::Message() << "policy=" << static_cast<int>(policy));

    // Build via construct() — this writes every per-member cache.
    sri::GenericStorage build_storage;
    Index original(std::ref(build_storage), 512, 4.0f, policy);
    construct(original, this->config_);

    auto results_original = this->RunSearch(original);
    const auto report_original = original.GetSizeReport();
    const auto bytes_original = sdsl::size_in_bytes(original);

    // Re-load into a fresh index against fresh storage. Constructor
    // args (block_size, storing_factor, policy) must match — the cache
    // key embeds them.
    sri::GenericStorage reload_storage;
    Index reloaded(std::ref(reload_storage), 512, 4.0f, policy);
    reloaded.load(this->config_);

    auto results_reloaded = this->RunSearch(reloaded);
    const auto report_reloaded = reloaded.GetSizeReport();
    const auto bytes_reloaded = sdsl::size_in_bytes(reloaded);

    // Search results must match the hand-coded expectations and the
    // pre-reload run.
    ASSERT_EQ(results_original.size(), this->patterns_.size());
    for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
      EXPECT_THAT(results_original[i],
                  testing::ElementsAreArray(this->patterns_[i].docs))
          << "original / " << this->patterns_[i].pattern;
      EXPECT_THAT(results_reloaded[i],
                  testing::ElementsAreArray(this->patterns_[i].docs))
          << "reloaded / " << this->patterns_[i].pattern;
    }

    // Size report: same fields in the same order with the same byte
    // counts. SDSL rank/select supports rebuild deterministically from
    // the underlying bitvector during PDLTreeCore::rebindRankSelect, so
    // no slack is needed here.
    ASSERT_EQ(report_original.size(), report_reloaded.size());
    for (std::size_t i = 0; i < report_original.size(); ++i) {
      EXPECT_EQ(report_original[i].name, report_reloaded[i].name);
      EXPECT_EQ(report_original[i].bytes, report_reloaded[i].bytes)
          << "field=" << report_original[i].name;
    }
    EXPECT_EQ(bytes_original, bytes_reloaded);
  }
}

}  // namespace
