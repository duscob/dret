//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/14/26.
//
// Cache-key collision audit, Task 45 of docs/pdl_indexes_tasks.md.
// Regression net for the gotcha documented in
// `generic_storage_key_collision.md`: when two indexes share a
// GenericStorage and the same logical key but different TItem types,
// IndexBaseWithExternalStorage::loadItemPtr must include the type
// hash in the storage key (not just the on-disk filename) — otherwise
// the second index loads via a null pointer and SEGVs inside its
// component's load(). The dret/index_base.h:loadItemPtr fix appended
// the SDSL type-hash to the in-memory storage key when
// add_type_hash=true; this test exercises that fix across PDL's
// (variant, storage_policy, get_doc, block_size, storing_factor)
// axis space.
//
// Acceptance: all cache entries are distinct in memory and on disk;
// no SEGV, no silent reuse.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"
#include "dret/pdl/storage_policy.h"

#include "base_test.h"

namespace {

using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;
using DefaultAlphabet = dret::Alphabet<>;
using DefaultCount = sri::RIndexCount<ExternalGenericStorage, DefaultAlphabet>;

// Plain-DA, the base PDL Plain instantiation.
using PlainDA = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>;
// RP-DA, for the variant axis.
using RPDA = dret::pdl::DocListIdxPDLRP<ExternalGenericStorage>;
// Plain-SLP, for the get-doc axis.
using PlainSLP = dret::pdl::DocListIdxPDLPlain<
    ExternalGenericStorage,
    DefaultAlphabet,
    DefaultCount,
    dret::pdl::PDLGetDocsSLP<ExternalGenericStorage, DefaultAlphabet::int_width>>;

class DocListResultVector : public std::vector<std::size_t> {
 public:
  void operator()(std::size_t doc) { this->push_back(doc); }
  void operator()() {
    std::sort(this->begin(), this->end());
    this->erase(std::unique(this->begin(), this->end()), this->end());
  }
};

class PDLCacheCollisionTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  // Smallish 5-doc fixture matching pdl_search_test.cpp's pattern
  // class coverage. The exact docs returned are not the load-bearing
  // assertion here; what matters is that every distinct-tuple PDL
  // index loaded into the same storage returns the same result as
  // the brute-force baseline AND that an earlier index's results do
  // not change after later indexes share the storage.
  const std::string data_ =
      "abracadabra\1barbara\1albatross\1amanda\1lasagna\1";

  const std::vector<std::string> patterns_ = {
      "a", "ab", "ar", "an", "abra", "lasagna", "xyz",
  };

  template <typename TIdx>
  std::vector<std::size_t> Search(TIdx& idx, const std::string& pattern) {
    DocListResultVector result;
    idx.Search(pattern, std::ref(result));
    result();
    return std::vector<std::size_t>(result.begin(), result.end());
  }
};

TEST_F(PDLCacheCollisionTest, AllAxesSharedStorageNoCollision) {
  // Brute-force baseline.
  sri::GenericStorage brute_storage;
  dret::DocListIdxBrute<ExternalGenericStorage> brute(std::ref(brute_storage));
  construct(brute, this->config_);
  std::vector<std::vector<std::size_t>> brute_expected;
  brute_expected.reserve(this->patterns_.size());
  for (const auto& p : this->patterns_) {
    brute_expected.push_back(this->Search(brute, p));
  }

  // Single shared storage across every PDL variant. Pre-fix
  // behavior: loading the second variant would SEGV inside its
  // component's load(). Post-fix: each (key + type-hash) is a
  // distinct map entry, so all variants coexist.
  sri::GenericStorage shared;
  const auto storage_size_at_start = shared.size();

  // Base configuration: Plain / OccurrenceWeighted / DA / 512 / 4.
  PlainDA base(std::ref(shared), 512, 4.0f,
               dret::pdl::StoragePolicy::OccurrenceWeighted);
  construct(base, this->config_);

  // Search and capture the base index's per-pattern output. We will
  // re-run this search after every other variant has loaded into the
  // same storage; the result must not change (no silent overwrite).
  std::vector<std::vector<std::size_t>> base_initial;
  base_initial.reserve(this->patterns_.size());
  for (const auto& p : this->patterns_) {
    base_initial.push_back(this->Search(base, p));
    EXPECT_EQ(base_initial.back(),
              brute_expected[base_initial.size() - 1])
        << "base / pattern " << p;
  }

  // --- variant axis: RP instead of Plain ---
  RPDA rp(std::ref(shared), 512, 4.0f,
          dret::pdl::StoragePolicy::OccurrenceWeighted);
  construct(rp, this->config_);
  for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
    EXPECT_EQ(this->Search(rp, this->patterns_[i]), brute_expected[i])
        << "variant axis (RP) / pattern " << this->patterns_[i];
  }

  // --- storage_policy axis: StoreAllInternal instead of OccurrenceWeighted ---
  PlainDA policy(std::ref(shared), 512, 4.0f,
                 dret::pdl::StoragePolicy::StoreAllInternal);
  construct(policy, this->config_);
  for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
    EXPECT_EQ(this->Search(policy, this->patterns_[i]), brute_expected[i])
        << "policy axis (StoreAllInternal) / pattern " << this->patterns_[i];
  }

  // --- block_size axis: 256 instead of 512 ---
  PlainDA bs(std::ref(shared), 256, 4.0f,
             dret::pdl::StoragePolicy::OccurrenceWeighted);
  construct(bs, this->config_);
  for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
    EXPECT_EQ(this->Search(bs, this->patterns_[i]), brute_expected[i])
        << "block_size axis (256) / pattern " << this->patterns_[i];
  }

  // --- storing_factor axis: 8 instead of 4 ---
  PlainDA sf(std::ref(shared), 512, 8.0f,
             dret::pdl::StoragePolicy::OccurrenceWeighted);
  construct(sf, this->config_);
  for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
    EXPECT_EQ(this->Search(sf, this->patterns_[i]), brute_expected[i])
        << "storing_factor axis (8) / pattern " << this->patterns_[i];
  }

  // --- get_doc axis: SLP instead of DA. PlainSLP's get-doc cache
  // entry has a different type than PlainDA's, so a class_to_hash
  // disambiguation is required at storage-key time even when both
  // indexes share the exact same logical key prefix. ---
  PlainSLP getdoc(std::ref(shared), 512, 4.0f,
                  dret::pdl::StoragePolicy::OccurrenceWeighted);
  construct(getdoc, this->config_);
  for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
    EXPECT_EQ(this->Search(getdoc, this->patterns_[i]), brute_expected[i])
        << "get_doc axis (SLP) / pattern " << this->patterns_[i];
  }

  // After every variant has loaded into `shared`, re-run the base
  // index's search. If any later variant's load had silently
  // overwritten the base's typed entry (the pre-fix collision class),
  // base would either SEGV or report wrong docs.
  for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
    EXPECT_EQ(this->Search(base, this->patterns_[i]), base_initial[i])
        << "base output changed after sharing storage; pattern "
        << this->patterns_[i];
  }

  // Distinct-in-memory check: every variant added at least one new
  // entry (often several — count_idx, get_docs, core under typed
  // keys). Six tuples total, so at least six entries beyond start.
  EXPECT_GE(shared.size() - storage_size_at_start, 6u)
      << "expected at least one new typed entry per tuple in shared storage";
}

}  // namespace
