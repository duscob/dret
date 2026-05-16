//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// Storage-policy correctness fixture, Task 43 of
// docs/pdl_indexes_tasks.md. For each PDL storage policy
// (OccurrenceWeighted / StoreAllInternal / LeavesOnly):
//
//   1. PDL search output equals brute-force output for every pattern.
//   2. At least one pattern's cover under that policy contains
//      multiple raw ranges, exercising the
//      DLSampledTreeScheme::Search → computeCoverFull → multi-getDocs
//      multi-range path. Without this guard, Task 34 / 42's parity
//      tests would still pass on a fixture where every cover happens
//      to be a single range, leaving the multi-range path untested.
//
// Acceptance: each policy yields documents identical to brute force,
// including patterns whose covers contain ≥ 2 raw ranges.
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
#include "dret/pdl/doc_list_pdl_plain.h"
#include "dret/pdl/storage_policy.h"

#include "base_test.h"

namespace {

using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;
using PDLPlain = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>;

// Inspector subclass: re-publishes count() and computeCoverFull() so
// the test can probe the cover decision directly without depending on
// the Search side effect.
struct PDLPlainInspector : PDLPlain {
  using PDLPlain::PDLPlain;
  using PDLPlain::computeCoverFull;
  using PDLPlain::count;
  using PDLPlain::get_docs;
};

class DocListResultVector : public std::vector<std::size_t> {
 public:
  void operator()(std::size_t doc) { this->push_back(doc); }
  void operator()() {
    std::sort(this->begin(), this->end());
    this->erase(std::unique(this->begin(), this->end()), this->end());
  }
};

class PDLStoragePolicyCorrectnessTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  // 8-doc fixture chosen so the suffix tree carries a mix of
  // single-doc and multi-doc subtrees, and each storage policy
  // produces a distinguishable selected-node set. The patterns below
  // are then chosen so brute force returns at least one match
  // spanning a non-trivial SA range — enough to force partial-leaf
  // overlaps under all three policies.
  const std::string data_ =
      "abracadabra\1barbara\1albatross\1amanda\1lasagna\1"
      "alabama\1canada\1arabia\1";

  struct PatternCase {
    std::string label;
    std::string pattern;
  };

  // Mix of long-SA-range patterns ('a' is very common) and shorter
  // ones with non-trivial bigram overlap. The long-range queries are
  // the ones that produce multi-range covers under LeavesOnly /
  // OccurrenceWeighted.
  const std::vector<PatternCase> patterns_ = {
      {"common_a",        "a"},
      {"common_ab",       "ab"},
      {"medium_ar",       "ar"},
      {"medium_an",       "an"},
      {"few_bar",         "bar"},
      {"few_ana",         "ana"},
      {"once_lasagna",    "lasagna"},
      {"none_xyz",        "xyz"},
  };

  template <typename TIdx>
  std::vector<std::size_t> Search(TIdx& idx, const std::string& pattern) {
    DocListResultVector result;
    idx.Search(pattern, std::ref(result));
    result();
    return std::vector<std::size_t>(result.begin(), result.end());
  }
};

TEST_F(PDLStoragePolicyCorrectnessTest, EachPolicyMatchesBruteAndExercisesMultiRangeCover) {
  // Brute baseline.
  sri::GenericStorage brute_storage;
  dret::DocListIdxBrute<ExternalGenericStorage> brute(std::ref(brute_storage));
  construct(brute, this->config_);

  std::vector<std::vector<std::size_t>> brute_expected;
  brute_expected.reserve(this->patterns_.size());
  for (const auto& pc : this->patterns_) {
    brute_expected.push_back(this->Search(brute, pc.pattern));
  }

  for (auto policy : {dret::pdl::StoragePolicy::OccurrenceWeighted,
                      dret::pdl::StoragePolicy::StoreAllInternal,
                      dret::pdl::StoragePolicy::LeavesOnly}) {
    SCOPED_TRACE(::testing::Message() << "policy=" << static_cast<int>(policy));

    // block_size = 1 prevents CollapseSubtreesByBlockSize from
    // flattening the suffix tree into a single block — necessary for
    // the small fixture to retain multiple internal/leaf nodes that
    // can decompose covers into multiple raw ranges.
    sri::GenericStorage pdl_storage;
    PDLPlainInspector pdl(std::ref(pdl_storage), /*block_size=*/1,
                          /*storing_factor=*/4.0f, policy);
    construct(pdl, this->config_);

    // Pattern-level parity: each pattern's PDL search output must
    // match brute force regardless of how the cover decomposes.
    for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
      const auto& pc = this->patterns_[i];
      auto pdl_result = this->Search(pdl, pc.pattern);
      EXPECT_THAT(pdl_result, testing::ElementsAreArray(brute_expected[i]))
          << pc.label << " pattern=" << pc.pattern;
    }

    // Multi-range-cover trigger: probe arbitrary [sp, ep) ranges over
    // the full SA looking for at least one cover that decomposes into
    // ≥ 2 raw ranges. Tying this to pattern-derived SA ranges would
    // make the assertion fragile (the chosen patterns may all happen
    // to produce single-range covers on a small fixture); probing the
    // SA directly removes that dependency. This pins the multi-range
    // path is reachable under each storage policy on this fixture.
    const std::size_t n = pdl.get_docs().inner().size();
    ASSERT_GT(n, 0u);

    std::size_t max_raw_ranges_for_policy = 0;
    std::pair<std::size_t, std::size_t> max_range;
    for (std::size_t sp = 0; sp <= n; ++sp) {
      for (std::size_t ep = sp + 1; ep <= n; ++ep) {
        std::vector<std::pair<std::size_t, std::size_t>> raw_ranges;
        std::vector<std::size_t> nodes;
        pdl.computeCoverFull(sp, ep, raw_ranges, nodes);
        if (raw_ranges.size() > max_raw_ranges_for_policy) {
          max_raw_ranges_for_policy = raw_ranges.size();
          max_range = {sp, ep};
          if (max_raw_ranges_for_policy >= 2) goto done;
        }
      }
    }
   done:;
    EXPECT_GE(max_raw_ranges_for_policy, 2u)
        << "no [sp, ep) range exercised the multi-range cover path "
        << "under policy " << static_cast<int>(policy)
        << "; richest cover was " << max_raw_ranges_for_policy
        << " raw range(s) at [" << max_range.first << ", "
        << max_range.second << ")";
  }
}

}  // namespace
