//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/12/26.
//
// DocListIdxPDLPlain smoke tests (Task 23 acceptance: class compiles,
// serializes/loads, and passes basic DA-backed search tests). Full
// typed parity vs brute force across all PDL variants lands in Tasks
// 33-34, which add the index to test/doc_list_test.cpp's existing
// type lists.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "dret/doc_list_index_brute.h"
#include "dret/pdl/doc_list_pdl_plain.h"

#include "base_test.h"

namespace {

using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;

class DocListIdxPDLPlainSmokeTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  sri::GenericStorage storage_;

  const std::string data_ = "TATA\1LATA\1LALA\1";

  struct SearchData {
    std::string pattern;
    std::vector<std::size_t> docs;
  };

  const std::vector<SearchData> search_data_ = {
      {"TAT", {0}},
      {"LAT", {1}},
      {"LAL", {2}},
      {"TA",  {0, 1}},
      {"LA",  {1, 2}},
      {"A",   {0, 1, 2}},
      {"TAL", {}},
  };
};

class DocListResultVector : public std::vector<std::size_t> {
 public:
  void operator()(std::size_t doc) { this->push_back(doc); }
  void operator()() {
    std::sort(this->begin(), this->end());
    this->erase(std::unique(this->begin(), this->end()), this->end());
  }
};

// Build a PDLPlain index over the fixture, search every pattern, compare
// to the hand-coded expected document lists.
TEST_F(DocListIdxPDLPlainSmokeTest, SearchMatchesExpected) {
  using Index = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>;
  Index index(std::ref(storage_));
  construct(index, this->config_);

  for (const auto& [pattern, docs] : this->search_data_) {
    DocListResultVector result;
    index.Search(pattern, std::ref(result));
    result();

    EXPECT_EQ(result.size(), docs.size()) << pattern;
    EXPECT_THAT(result, testing::ElementsAreArray(docs)) << pattern;
  }
}

// Task 23 "serializes, loads" acceptance: construct once (which writes
// the PDLTreeCore + count_idx + DA caches), then build a fresh Index
// instance against an empty GenericStorage and call load() — that
// exercises the disk round-trip via PDLTreeCore::load,
// RIndexCount::load, and GetDocDA::load. Cannot use
// sdsl::store_to_cache(index, ...) here because the index's type-hash
// path constructs Index{}, which is ill-formed when TStorage is
// std::reference_wrapper.
TEST_F(DocListIdxPDLPlainSmokeTest, SerializeLoadRoundTripPreservesSearch) {
  using Index = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>;

  // First instance: drives construct() to write all per-member caches.
  Index original(std::ref(storage_));
  construct(original, this->config_);

  // Reload into a fresh instance against fresh storage.
  sri::GenericStorage fresh_storage;
  Index reloaded(std::ref(fresh_storage));
  reloaded.load(this->config_);

  for (const auto& [pattern, docs] : this->search_data_) {
    DocListResultVector result;
    reloaded.Search(pattern, std::ref(result));
    result();

    EXPECT_EQ(result.size(), docs.size()) << pattern;
    EXPECT_THAT(result, testing::ElementsAreArray(docs)) << pattern;
  }
}

// Compare PDLPlain against the brute-force baseline on the same fixture
// for an extra safety net beyond hand-coded expectations.
TEST_F(DocListIdxPDLPlainSmokeTest, MatchesBruteForce) {
  using PDLIndex = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>;
  using BruteIndex = dret::DocListIdxBrute<ExternalGenericStorage>;

  PDLIndex pdl(std::ref(storage_));
  construct(pdl, this->config_);

  // Brute force uses its own storage to keep independent caches.
  sri::GenericStorage brute_storage;
  BruteIndex brute(std::ref(brute_storage));
  construct(brute, this->config_);

  for (const auto& [pattern, _] : this->search_data_) {
    DocListResultVector pdl_result;
    pdl.Search(pattern, std::ref(pdl_result));
    pdl_result();

    DocListResultVector brute_result;
    brute.Search(pattern, std::ref(brute_result));
    brute_result();

    EXPECT_THAT(pdl_result, testing::ElementsAreArray(brute_result)) << pattern;
  }
}

}  // namespace
