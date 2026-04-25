//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/vlc_vector.hpp>

#include "dret/basic_slp_span_length.h"
#include "dret/differential_slp.h"
#include "dret/doc_list_index_brute.h"
#include "dret/doc_list_index_rmq.h"
#include "dret/doc_list_sampled_tree_dgcda.h"
#include "dret/doc_list_sampled_tree_gcda.h"

#include "base_test.h"

//~~~~~~~


template <typename TIndex>
class DocListIndexConstructTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  const std::string data_ = "MINIMUM\1MINIMAL\1MINIMIZES\1";
};

using DocListIndexConstructTypes = ::testing::Types<
    dret::DocListIdxBrute<>,
    dret::gcda::DocListIdxGCDA<>,
    dret::dgcda::DocListIdxDGCDA<>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                                 dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                                 dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                             dret::rmq::SadaCore<dret::GenericStorage>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                             dret::rmq::IlcpCore<dret::GenericStorage>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                             dret::rmq::CilcpCore<dret::GenericStorage>>>;

TYPED_TEST_SUITE(DocListIndexConstructTypedTests, DocListIndexConstructTypes);

TYPED_TEST(DocListIndexConstructTypedTests, construct) {
  auto key_index = "index";
  {
    TypeParam index;
    construct(index, this->config_);
    sdsl::store_to_cache(index, key_index, this->config_, true);
  }

  TypeParam index;
  sdsl::load_from_cache(index, key_index, this->config_, true);
}

//~~~~~~~


using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;

template <typename TIndex>
class DocListIndexSearchTypedTests : public BaseConfigTests<8> {
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
      {"TA", {0, 1}},
      {"LA", {1, 2}},
      {"A", {0, 1, 2}},
      {"TAL", {}},
  };
};

using DocListIndexSearchTypes = ::testing::Types<
    dret::DocListIdxBrute<ExternalGenericStorage>,
    dret::gcda::DocListIdxGCDA<ExternalGenericStorage>,
    dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage>,
    dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage,
                                 dret::Alphabet<>,
                                 sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                                 dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
                                 dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::SrIdxGeneric<sri::SrIndexValidArea<ExternalGenericStorage, dret::Alphabet<>>, 16>,
                             dret::rmq::SadaCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::SrIdxGeneric<sri::SrIndexValidArea<ExternalGenericStorage, dret::Alphabet<>>, 16>,
                             dret::rmq::IlcpCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::SrIdxGeneric<sri::SrIndexValidArea<ExternalGenericStorage, dret::Alphabet<>>, 16>,
                             dret::rmq::CilcpCore<ExternalGenericStorage>>>;

TYPED_TEST_SUITE(DocListIndexSearchTypedTests, DocListIndexSearchTypes);

class DocListResult {
 public:
  virtual ~DocListResult() = default;

  virtual void operator()(std::size_t t_doc) = 0;

  virtual void operator()() {}

  virtual void Print(std::ostream& t_os) const = 0;
};

class DocListResultVector : public DocListResult, public std::vector<std::size_t> {
 public:
  void operator()(std::size_t t_doc) override {
    emplace_back(t_doc);
  }

  using std::vector<std::size_t>::begin;
  using std::vector<std::size_t>::end;

  void operator()() override {
    sort(begin(), end());
    erase(unique(begin(), end()), end());
  }

  void Print(std::ostream& t_os) const override {
    for (const auto& item : *this) {
      t_os << item << '\n';
    }
  }
};

TYPED_TEST(DocListIndexSearchTypedTests, search) {
  TypeParam index(std::ref(this->storage_));
  construct(index, this->config_);

  for (const auto& [pattern, docs] : this->search_data_) {
    DocListResultVector result;
    index.Search(pattern, std::ref(result));
    result();

    EXPECT_EQ(result.size(), docs.size()) << pattern;
    EXPECT_THAT(result, testing::ElementsAreArray(docs)) << pattern;
  }
}

//~~~~~~~


template <typename TDSLP>
class DifferentialSLPExpandTypedTests : public ::testing::Test {};

using DifferentialSLPExpandTypes =
    ::testing::Types<dret::DifferentialSLP<>,
                     dret::DifferentialSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<sdsl::int_vector<>>>>,
                     dret::DifferentialSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<sdsl::int_vector<>>>>,
                     dret::DifferentialSLP<grammar::SLP<>,
                                           sdsl::enc_vector<>,
                                           sdsl::enc_vector<>,
                                           sdsl::enc_vector<>>,
                     dret::DifferentialSLP<grammar::SLP<>,
                                           sdsl::dac_vector<>,
                                           sdsl::dac_vector<>,
                                           sdsl::dac_vector<>>,
                     dret::DifferentialSLP<grammar::SLP<>,
                                           sdsl::vlc_vector<>,
                                           sdsl::vlc_vector<>,
                                           sdsl::vlc_vector<>>>;

TYPED_TEST_SUITE(DifferentialSLPExpandTypedTests, DifferentialSLPExpandTypes);

TYPED_TEST(DifferentialSLPExpandTypedTests, expand_roundtrips_da) {
  const std::vector<std::uint64_t> da_vec = {
      0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 3, 3, 3, 3, 3, 1, 2, 0, 1, 2, 0, 1, 2,
  };
  sdsl::int_vector<> da(da_vec.size(), 0, 64);
  for (std::size_t i = 0; i < da_vec.size(); ++i)
    da[i] = da_vec[i];
  sdsl::util::bit_compress(da);

  TypeParam dslp;
  dslp.Compute(da, /*block_size=*/4);

  std::vector<std::size_t> expanded;
  auto report = [&expanded](auto value) {
    expanded.emplace_back(static_cast<std::size_t>(value));
  };
  dret::ExpandSLP(dslp, 0, da.size(), report);

  ASSERT_EQ(expanded.size(), da_vec.size());
  for (std::size_t i = 0; i < da_vec.size(); ++i)
    EXPECT_EQ(expanded[i], da_vec[i]) << "position " << i;
}
