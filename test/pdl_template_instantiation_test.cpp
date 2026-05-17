//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// Index-level template-instantiation tests for the bitvector axis,
// Task 37 of docs/pdl_indexes_tasks.md. Task 18's pdl_template_test.cpp
// already covers PDLTreeCore<TBitvector> and BCCodec<TBitvector> in
// isolation with sdsl::sd_vector<> and sdsl::bit_vector. This file
// extends that coverage to the full DocListIdxPDL{Plain,RP,BC} index
// templates: each PDL index gets two typed instantiations (default
// sd_vector tree-core marker + an explicit bit_vector tree-core
// marker), and both must compile, build via construct(), and return
// correct search results.
//
// Acceptance: minimal runtime tests pass for both bitvector families.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/sd_vector.hpp>

#include "dret/config.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"
#include "dret/pdl/set_codecs.h"
#include "sr-index/r_index.h"

#include "base_test.h"

namespace {

using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;
using DefaultAlphabet = dret::Alphabet<>;
using DefaultCount = sri::RIndexCount<ExternalGenericStorage, DefaultAlphabet>;
using DefaultGetDocs = dret::pdl::PDLGetDocsDA<ExternalGenericStorage,
                                               DefaultAlphabet::int_width>;

class DocListResultVector : public std::vector<std::size_t> {
 public:
  void operator()(std::size_t doc) { this->push_back(doc); }
  void operator()() {
    std::sort(this->begin(), this->end());
    this->erase(std::unique(this->begin(), this->end()), this->end());
  }
};

// One typed-test traits struct per (codec, bitvector) instantiation.
template <typename TIdx>
struct InstTraits {
  using Index = TIdx;
};

// The alias templates DocListIdxPDL{Plain,RP,BC} lock the codec at the 5th
// position of the underlying DocListIdxPDL template, so each alias takes 7
// parameters — the codec is implicit in the alias name.

// --- Plain ---
using PlainSdVector =
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage,
                                  DefaultAlphabet,
                                  DefaultCount,
                                  DefaultGetDocs,
                                  sdsl::sd_vector<>>;
using PlainBitVector =
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage,
                                  DefaultAlphabet,
                                  DefaultCount,
                                  DefaultGetDocs,
                                  sdsl::bit_vector>;

// --- RP ---
using RPSdVector =
    dret::pdl::DocListIdxPDLRP<ExternalGenericStorage,
                               DefaultAlphabet,
                               DefaultCount,
                               DefaultGetDocs,
                               sdsl::sd_vector<>>;
using RPBitVector =
    dret::pdl::DocListIdxPDLRP<ExternalGenericStorage,
                               DefaultAlphabet,
                               DefaultCount,
                               DefaultGetDocs,
                               sdsl::bit_vector>;

// --- BC ---
using BCSdVector =
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage,
                               DefaultAlphabet,
                               DefaultCount,
                               DefaultGetDocs,
                               sdsl::sd_vector<>>;
using BCBitVector =
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage,
                               DefaultAlphabet,
                               DefaultCount,
                               DefaultGetDocs,
                               sdsl::bit_vector>;

template <typename TTraits>
class PDLIndexBitvectorTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  sri::GenericStorage storage_;
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
};

using PDLIndexBitvectorTraitsList = ::testing::Types<
    InstTraits<PlainSdVector>,
    InstTraits<PlainBitVector>,
    InstTraits<RPSdVector>,
    InstTraits<RPBitVector>,
    InstTraits<BCSdVector>,
    InstTraits<BCBitVector>>;

TYPED_TEST_SUITE(PDLIndexBitvectorTypedTests, PDLIndexBitvectorTraitsList);

// Compile + smoke: construct the index against the small fixture and
// run the canonical pattern set. If the bitvector swap broke the
// PDLTreeCore tree-core marker or its rank/select rebind, this
// surfaces immediately as wrong results (or a build-time failure).
TYPED_TEST(PDLIndexBitvectorTypedTests, CompilesAndSearchesCorrectly) {
  using Index = typename TypeParam::Index;
  Index index(std::ref(this->storage_));
  construct(index, this->config_);

  for (const auto& [pattern, expected] : this->patterns_) {
    SCOPED_TRACE(::testing::Message() << "pattern=" << pattern);
    DocListResultVector result;
    index.Search(pattern, std::ref(result));
    result();

    EXPECT_THAT(result, testing::ElementsAreArray(expected));
  }
}

}  // namespace
