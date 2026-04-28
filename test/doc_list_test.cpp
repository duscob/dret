//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/vlc_vector.hpp>

#include "dret/basic_slp_span_length.h"
#include "dret/differential_slp.h"
#include "dret/doc_list_index_brute.h"
#include "dret/doc_list_index_rmq.h"
#include "dret/doc_list_rmq_scheme.h"
#include "dret/doc_list_sampled_tree_dgcda.h"
#include "dret/doc_list_idx_slp.h"
#include "dret/doc_list_sampled_tree_gcda.h"

#include "base_test.h"

//~~~~~~~

template <typename TStorage>
using RMQCountIdx = sri::RIndexCount<TStorage, dret::Alphabet<>>;

template <typename TStorage>
using RMQGetDocSLP = dret::rmq::GetDocSLP<TStorage>;

template <typename TStorage>
using RMQGetDocDSLP = dret::rmq::GetDocDSLP<TStorage>;

template <typename TStorage>
using RMQSadaSLPCore = dret::rmq::SadaCore<TStorage,
                                           dret::Alphabet<>::int_width,
                                           sdsl::rmq_succinct_sct<true>,
                                           sdsl::sd_vector<>,
                                           RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQSadaDSLPCore = dret::rmq::SadaCore<TStorage,
                                            dret::Alphabet<>::int_width,
                                            sdsl::rmq_succinct_sct<true>,
                                            sdsl::sd_vector<>,
                                            RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPCore = dret::rmq::IlcpCore<TStorage,
                                           dret::Alphabet<>::int_width,
                                           sdsl::bit_vector,
                                           sdsl::rmq_succinct_sct<true>,
                                           sdsl::sd_vector<>,
                                           RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQIlcpDSLPCore = dret::rmq::IlcpCore<TStorage,
                                            dret::Alphabet<>::int_width,
                                            sdsl::bit_vector,
                                            sdsl::rmq_succinct_sct<true>,
                                            sdsl::sd_vector<>,
                                            RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPCore = dret::rmq::CilcpCore<TStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::bit_vector,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQCilcpDSLPCore = dret::rmq::CilcpCore<TStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::bit_vector,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQSadaSLPIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQSadaSLPCore<TStorage>>;

template <typename TStorage>
using RMQSadaDSLPIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQSadaDSLPCore<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQIlcpSLPCore<TStorage>>;

template <typename TStorage>
using RMQIlcpDSLPIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQIlcpDSLPCore<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQCilcpSLPCore<TStorage>>;

template <typename TStorage>
using RMQCilcpDSLPIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQCilcpDSLPCore<TStorage>>;

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
    dret::gcda::DocListIdxGCDA<dret::GenericStorage,
                               dret::Alphabet<>,
                               sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                               grammar::CompactBPSLP<>>,
    dret::gcda::DocListIdxGCDA<dret::GenericStorage,
                               dret::Alphabet<>,
                               sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                               grammar::CompactLOUDSSLP<>>,
    dret::gcda::DocListIdxGCDA<dret::GenericStorage,
                               dret::Alphabet<>,
                               sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                               grammar::CombinedSLPWithUnitCover<>>,
    dret::DocListIdxSLP<>,
    dret::dgcda::DocListIdxDGCDA<>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaCore<dret::GenericStorage>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpCore<dret::GenericStorage>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpCore<dret::GenericStorage>>,
    RMQSadaSLPIndex<dret::GenericStorage>,
    RMQIlcpSLPIndex<dret::GenericStorage>,
    RMQCilcpSLPIndex<dret::GenericStorage>,
    RMQSadaDSLPIndex<dret::GenericStorage>,
    RMQIlcpDSLPIndex<dret::GenericStorage>,
    RMQCilcpDSLPIndex<dret::GenericStorage>>;

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
    dret::gcda::DocListIdxGCDA<ExternalGenericStorage,
                               dret::Alphabet<>,
                               sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                               grammar::CompactBPSLP<>>,
    dret::gcda::DocListIdxGCDA<ExternalGenericStorage,
                               dret::Alphabet<>,
                               sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                               grammar::CompactLOUDSSLP<>>,
    dret::gcda::DocListIdxGCDA<ExternalGenericStorage,
                               dret::Alphabet<>,
                               sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                               grammar::CombinedSLPWithUnitCover<>>,
    dret::DocListIdxSLP<ExternalGenericStorage>,
    dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage>,
    dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpCore<ExternalGenericStorage>>,
    RMQSadaSLPIndex<ExternalGenericStorage>,
    RMQIlcpSLPIndex<ExternalGenericStorage>,
    RMQCilcpSLPIndex<ExternalGenericStorage>,
    RMQSadaDSLPIndex<ExternalGenericStorage>,
    RMQIlcpDSLPIndex<ExternalGenericStorage>,
    RMQCilcpDSLPIndex<ExternalGenericStorage>>;

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

class RMQSLPCacheReuseTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  sri::GenericStorage storage_;
  const std::string data_ = "TATA\1LATA\1LALA\1";
};

TEST_F(RMQSLPCacheReuseTest, rmq_slp_reuses_gcda_slp_cache) {
  using GetDocSLP = RMQGetDocSLP<ExternalGenericStorage>;
  using TSLP = typename GetDocSLP::SLP;
  using RMQSadaSLP = RMQSadaSLPIndex<ExternalGenericStorage>;

  dret::gcda::DocListIdxGCDA<ExternalGenericStorage> gcda(std::ref(storage_), 512, 4);
  construct(gcda, config_);

  const auto key_slp = std::format("512-4_{}", config_.keys[dret::conf::kGCDA][dret::conf::kSLP].get<std::string>());
  const auto slp_path = sdsl::cache_file_name<TSLP>(key_slp, config_);
  ASSERT_TRUE(std::filesystem::exists(slp_path));
  const auto before_time = std::filesystem::last_write_time(slp_path);
  const auto before_size = std::filesystem::file_size(slp_path);

  RMQSadaSLP rmq_sada_slp(std::ref(storage_));
  construct(rmq_sada_slp, config_);

  ASSERT_TRUE(std::filesystem::exists(slp_path));
  EXPECT_EQ(std::filesystem::file_size(slp_path), before_size);
  EXPECT_EQ(std::filesystem::last_write_time(slp_path), before_time);
}

TEST_F(RMQSLPCacheReuseTest, rmq_dslp_reuses_dgcda_dslp_cache) {
  using GetDocDSLP = RMQGetDocDSLP<ExternalGenericStorage>;
  using TDSLP = typename GetDocDSLP::DSLP;
  using RMQSadaDSLP = RMQSadaDSLPIndex<ExternalGenericStorage>;

  dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage> dgcda(std::ref(storage_), 512, 4);
  construct(dgcda, config_);

  const auto key_dslp = std::format("512-4_{}", config_.keys[dret::conf::kDGCDA][dret::conf::kSLP].get<std::string>());
  const auto dslp_path = sdsl::cache_file_name<TDSLP>(key_dslp, config_);
  ASSERT_TRUE(std::filesystem::exists(dslp_path));
  const auto before_time = std::filesystem::last_write_time(dslp_path);
  const auto before_size = std::filesystem::file_size(dslp_path);

  RMQSadaDSLP rmq_sada_dslp(std::ref(storage_));
  construct(rmq_sada_dslp, config_);

  ASSERT_TRUE(std::filesystem::exists(dslp_path));
  EXPECT_EQ(std::filesystem::file_size(dslp_path), before_size);
  EXPECT_EQ(std::filesystem::last_write_time(dslp_path), before_time);
}

//~~~~~~~


template <typename TDSLP>
class DifferentialSLPExpandTypedTests : public ::testing::Test {};

using DifferentialSLPExpandTypes =
    ::testing::Types<dret::DifferentialSLP<>,
                     dret::DifferentialSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<sdsl::int_vector<>>>>,
                     dret::DifferentialSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<sdsl::int_vector<>>>>,
                     dret::DifferentialSLP<grammar::SLP<>, sdsl::enc_vector<>, sdsl::enc_vector<>, sdsl::enc_vector<>>,
                     dret::DifferentialSLP<grammar::SLP<>, sdsl::dac_vector<>, sdsl::dac_vector<>, sdsl::dac_vector<>>,
                     dret::DifferentialSLP<grammar::SLP<>, sdsl::vlc_vector<>, sdsl::vlc_vector<>, sdsl::vlc_vector<>>>;

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

//~~~~~~~
// ListDocsRMQScheme unit tests
//
// Uses a synthetic DA and the SADA prev_doc array to exercise the stateless
// scheme directly — no full index construction required.


class ListDocsRMQSchemeTest : public ::testing::Test {
 protected:
  // DA = {0, 1, 2, 0, 1}  (3 distinct docs, 5 SA positions)
  // prev_doc[i] = last_occ[da[i]] before position i (0 when first occurrence).
  // prev_doc = {0, 0, 0, 0, 1}
  void SetUp() override {
    const std::vector<std::size_t> da_vec = {0, 1, 2, 0, 1};
    n_doc_ = 3;

    sdsl::int_vector<> prev_doc(da_vec.size(), 0, sdsl::bits::hi(da_vec.size()) + 1);
    std::vector<std::size_t> last_occ(n_doc_, 0);
    for (std::size_t i = 0; i < da_vec.size(); ++i) {
      prev_doc[i] = last_occ[da_vec[i]];
      last_occ[da_vec[i]] = i;
    }

    da_.resize(da_vec.size());
    for (std::size_t i = 0; i < da_vec.size(); ++i)
      da_[i] = da_vec[i];

    rmq_ = sdsl::rmq_succinct_sct<true>(&prev_doc);
  }

  std::vector<std::size_t> query(std::size_t sp, std::size_t ep_closed) {
    dret::rmq::MarkedReported mr(n_doc_);
    std::vector<std::size_t> reported;
    auto get_doc = [this](std::size_t k) {
      return static_cast<std::size_t>(da_[k]);
    };
    auto report = [&mr, &reported](std::size_t /*k*/, std::size_t d) {
      mr.mark(d);
      reported.push_back(d);
    };
    dret::rmq::ListDocsRMQScheme(sp, ep_closed + 1, rmq_, get_doc, mr, report);
    std::sort(reported.begin(), reported.end());
    return reported;
  }

  sdsl::int_vector<> da_;
  sdsl::rmq_succinct_sct<true> rmq_;
  std::size_t n_doc_ = 0;
};

TEST_F(ListDocsRMQSchemeTest, full_range_reports_all_docs) {
  EXPECT_EQ(query(0, 4), (std::vector<std::size_t>{0, 1, 2}));
}

TEST_F(ListDocsRMQSchemeTest, suffix_range_reports_two_docs) {
  EXPECT_EQ(query(3, 4), (std::vector<std::size_t>{0, 1}));
}

TEST_F(ListDocsRMQSchemeTest, single_position_reports_one_doc) {
  EXPECT_EQ(query(4, 4), (std::vector<std::size_t>{1}));
}

TEST_F(ListDocsRMQSchemeTest, middle_range_deduplicates) {
  // DA[1]=1, DA[2]=2, DA[3]=0
  EXPECT_EQ(query(1, 3), (std::vector<std::size_t>{0, 1, 2}));
}

TEST_F(ListDocsRMQSchemeTest, empty_range_reports_nothing) {
  EXPECT_EQ(query(2, 1), (std::vector<std::size_t>{}));
}
