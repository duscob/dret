//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// Search-correctness tests for every PDL variant (Plain / RP / BC)
// crossed with every storage policy
// (OccurrenceWeighted / StoreAllInternal / LeavesOnly), Task 34 of
// docs/pdl_indexes_tasks.md.
//
// Each typed case builds the PDL index AND a brute-force baseline on
// the same fixture, then asserts identical sorted/deduplicated zero-
// based doc-id output for patterns that occur (a) exactly once,
// (b) many times spread across multiple docs (with intra-doc
// repetition), (c) in every document, and (d) nowhere. This catches
// regressions in both the tree-cover path and the stored-set codec
// expansion path that the round-trip test in
// pdl_construct_load_test.cpp would not surface.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "dret/doc_list_index_brute.h"
#include "dret/pdl/doc_list_pdl_bc.h"
#include "dret/pdl/doc_list_pdl_plain.h"
#include "dret/pdl/doc_list_pdl_rp.h"
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
class DocListIdxPDLSearchTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  // Five docs chosen so the `pattern_kinds_` table below covers each
  // occurrence class:
  //   doc 0: "abracadabra"   — 'a' x5, 'abra' x2, unique 'cad'
  //   doc 1: "barbara"       — 'a' x3, 'bar' x2
  //   doc 2: "albatross"     — 'a' x2, 'al' x1, unique 'tross'
  //   doc 3: "amanda"        — 'a' x3, 'an' x1
  //   doc 4: "lasagna"       — 'a' x3, 'las' x1
  // Document delimiter is byte 0x01 (matches the 8-bit fixture in
  // doc_list_test.cpp / pdl_construct_load_test.cpp).
  const std::string data_ =
      "abracadabra\1barbara\1albatross\1amanda\1lasagna\1";

  // Patterns labelled by the occurrence class they exercise. The
  // "expected_docs" field is for documentation only — the actual
  // assertion compares PDL output to brute-force output, which is the
  // ground truth.
  struct PatternCase {
    std::string label;
    std::string pattern;
  };

  const std::vector<PatternCase> pattern_kinds_ = {
      // --- in every doc (5/5) ---
      {"in_all_letter_a",    "a"},
      // --- in many docs (2-4/5) ---
      {"many_substring_ar",  "ar"},
      {"many_substring_an",  "an"},
      // --- in exactly one doc ---
      {"once_unique_cad",    "cad"},
      {"once_unique_tross",  "tross"},
      {"once_full_doc_amanda", "amanda"},
      // --- not at all ---
      {"none_zzz",           "zzz"},
      {"none_xy",            "xy"},
  };

  template <typename TIdx>
  std::vector<std::size_t> Search(TIdx& idx, const std::string& pattern) {
    DocListResultVector result;
    idx.Search(pattern, std::ref(result));
    result();
    return std::vector<std::size_t>(result.begin(), result.end());
  }
};

using PDLSearchTypes = ::testing::Types<
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLRP<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage>>;

TYPED_TEST_SUITE(DocListIdxPDLSearchTypedTests, PDLSearchTypes);

TYPED_TEST(DocListIdxPDLSearchTypedTests, MatchesBruteForceAcrossAllPolicies) {
  using PDLIndex = TypeParam;
  using BruteIndex = dret::DocListIdxBrute<ExternalGenericStorage>;

  // Brute force baseline (independent of policy / codec).
  sri::GenericStorage brute_storage;
  BruteIndex brute(std::ref(brute_storage));
  construct(brute, this->config_);

  // Pre-compute brute-force expectations once per pattern.
  std::vector<std::vector<std::size_t>> brute_expected;
  brute_expected.reserve(this->pattern_kinds_.size());
  for (const auto& pc : this->pattern_kinds_) {
    brute_expected.push_back(this->Search(brute, pc.pattern));
  }

  for (auto policy : {dret::pdl::StoragePolicy::OccurrenceWeighted,
                      dret::pdl::StoragePolicy::StoreAllInternal,
                      dret::pdl::StoragePolicy::LeavesOnly}) {
    SCOPED_TRACE(::testing::Message() << "policy=" << static_cast<int>(policy));

    sri::GenericStorage pdl_storage;
    PDLIndex pdl(std::ref(pdl_storage), 512, 4.0f, policy);
    construct(pdl, this->config_);

    for (std::size_t i = 0; i < this->pattern_kinds_.size(); ++i) {
      const auto& pc = this->pattern_kinds_[i];
      auto pdl_result = this->Search(pdl, pc.pattern);
      EXPECT_THAT(pdl_result, testing::ElementsAreArray(brute_expected[i]))
          << pc.label << " pattern=" << pc.pattern;
    }
  }
}

}  // namespace
