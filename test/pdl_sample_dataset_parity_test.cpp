//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// PDL vs brute-force parity on a richer sample dataset, Task 42 of
// docs/pdl_indexes_tasks.md. Task 34's pdl_search_test.cpp validates
// each variant on a five-doc unit fixture; this test moves to a
// twenty-doc English-text fixture so the search path is exercised
// with realistic suffix-tree topology and a wider mix of pattern
// frequencies.
//
// The fixture is curated so that:
//   - Every doc contains the letter 'e', giving a true all-doc query
//     (and exercises the all-doc sentinel at codec build time).
//   - Several patterns occur in many but not all docs (most-common
//     case in real data).
//   - Several patterns occur in exactly one doc.
//   - Two patterns are intentionally absent (no-doc result).
//
// Acceptance: PDL variants match brute-force document sets exactly,
// with zero-based ids and no sentinel leakage.
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
class PDLSampleDatasetParityTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  // 20 English-text docs separated by 0x01. Every doc contains at
  // least one 'e', which makes the single-character pattern "e" an
  // all-doc query (5 occurrences in some docs, 1 in others — the
  // all-doc sentinel triggers at PDL build time).
  // Sources: classic short prose / nursery rhymes, varied lengths.
  const std::string data_ =
      "the quick brown fox jumps over the lazy dog\1"
      "to be or not to be that is the question\1"
      "all the worlds a stage and all the men and women merely players\1"
      "she sells seashells by the seashore\1"
      "peter piper picked a peck of pickled peppers\1"
      "humpty dumpty sat on a wall and had a great fall\1"
      "the rain in spain stays mainly in the plain\1"
      "how much wood could a woodchuck chuck before he had to stop\1"
      "i never saw a purple cow i never hope to see one\1"
      "fuzzy wuzzy was a bear who had no hair he was hardly fuzzy\1"
      "betty bought a bit of butter the bitter butter made her mad\1"
      "a stitch in time saves nine the early bird catches the worm\1"
      "row row row your boat gently down the stream\1"
      "twinkle twinkle little star how i wonder what you are\1"
      "mary had a little lamb whose fleece was white as snow\1"
      "old macdonald had a farm with a moo here and a quack there\1"
      "jack and jill went up the hill to fetch a pail of water\1"
      "the itsy bitsy spider crawled up the water spout\1"
      "hickory dickory dock the mouse ran up the clock\1"
      "london bridge is falling down my fair lady\1";

  struct PatternCase {
    std::string label;
    std::string pattern;
  };

  // Patterns curated so brute force's expected output covers every
  // occurrence class. The actual reference is brute-force output
  // captured at runtime; the labels are documentation only.
  const std::vector<PatternCase> patterns_ = {
      // --- in every doc (n_doc-wide; exercises all-doc sentinel) ---
      {"all_letter_e",        "e"},
      // --- in many docs (most-common shape in real data) ---
      {"many_the",            "the"},
      {"many_a_word",         " a "},
      {"many_letter_o",       "o"},
      // --- in a small handful of docs ---
      {"few_had",             "had"},
      {"few_water",           "water"},
      {"few_woodchuck",       "woodchuck"},
      // --- in exactly one doc ---
      {"once_seashells",      "seashells"},
      {"once_humpty",         "humpty"},
      {"once_stitch_in_time", "a stitch in time"},
      // --- not at any doc ---
      {"none_xyzzy",          "xyzzy"},
      {"none_kangaroo",       "kangaroo"},
  };

  template <typename TIdx>
  std::vector<std::size_t> Search(TIdx& idx, const std::string& pattern) {
    DocListResultVector result;
    idx.Search(pattern, std::ref(result));
    result();
    return std::vector<std::size_t>(result.begin(), result.end());
  }
};

using PDLSampleDatasetTypes = ::testing::Types<
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLRP<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage>>;

TYPED_TEST_SUITE(PDLSampleDatasetParityTypedTests, PDLSampleDatasetTypes);

TYPED_TEST(PDLSampleDatasetParityTypedTests, MatchesBruteForceAcrossPolicies) {
  using PDLIndex = TypeParam;
  using BruteIndex = dret::DocListIdxBrute<ExternalGenericStorage>;

  // Brute force baseline: build once, capture results per pattern.
  sri::GenericStorage brute_storage;
  BruteIndex brute(std::ref(brute_storage));
  construct(brute, this->config_);

  std::vector<std::vector<std::size_t>> brute_expected;
  brute_expected.reserve(this->patterns_.size());
  for (const auto& pc : this->patterns_) {
    brute_expected.push_back(this->Search(brute, pc.pattern));
  }

  // n_doc baseline: total doc count derived from the fixture itself
  // (count document-delimiter bytes 0x01). The codec sentinel is at
  // value n_doc, so legal doc ids must be strictly less.
  const std::size_t n_doc =
      std::count(this->data_.begin(), this->data_.end(), '\x01');
  ASSERT_GT(n_doc, 0u);

  for (auto policy : {dret::pdl::StoragePolicy::OccurrenceWeighted,
                      dret::pdl::StoragePolicy::StoreAllInternal,
                      dret::pdl::StoragePolicy::LeavesOnly}) {
    SCOPED_TRACE(::testing::Message() << "policy=" << static_cast<int>(policy));

    sri::GenericStorage pdl_storage;
    PDLIndex pdl(std::ref(pdl_storage), 512, 4.0f, policy);
    construct(pdl, this->config_);

    for (std::size_t i = 0; i < this->patterns_.size(); ++i) {
      const auto& pc = this->patterns_[i];
      auto pdl_result = this->Search(pdl, pc.pattern);

      // Exact match against brute force.
      EXPECT_THAT(pdl_result, testing::ElementsAreArray(brute_expected[i]))
          << pc.label << " pattern=\"" << pc.pattern << "\"";

      // Zero-based ids, no sentinel leakage. Every doc id must be
      // strictly less than n_doc, which catches the sentinel-leak
      // class of bug end-to-end (codec sentinel = n_doc).
      for (auto d : pdl_result) {
        EXPECT_LT(d, n_doc) << "sentinel leak in " << pc.label;
      }
    }
  }
}

}  // namespace
