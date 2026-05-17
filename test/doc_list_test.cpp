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

#include "dret/slp/basic_slp_span_length.h"
#include "dret/slp/differential_slp.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/rmq/doc_list_rmq_scheme.h"
#include "dret/doc_list/doc_list_slp.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"

#include "base_test.h"

//~~~~~~~

template <typename TStorage>
using RMQCountIdx = sri::RIndexCount<TStorage, dret::Alphabet<>>;

template <typename TStorage>
using RMQGetDocSLP = dret::rmq::GetDocSLP<TStorage>;

template <typename TStorage>
using RMQGetDocSLP_NS = dret::rmq::GetDocSLP_NS<TStorage>;

using BareSLP_Raw = grammar::SLP<>;
using BareSLP_DV = grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>;
using BareSLP_VV = grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>;

template <typename TStorage>
using RMQGetDocSLP_NS_Raw = dret::rmq::GetDocSLP_NS<TStorage, dret::Alphabet<>::int_width, BareSLP_Raw>;

template <typename TStorage>
using RMQGetDocSLP_NS_DV = dret::rmq::GetDocSLP_NS<TStorage, dret::Alphabet<>::int_width, BareSLP_DV>;

template <typename TStorage>
using RMQGetDocSLP_NS_VV = dret::rmq::GetDocSLP_NS<TStorage, dret::Alphabet<>::int_width, BareSLP_VV>;

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
                                           sdsl::sd_vector<>,
                                           sdsl::rmq_succinct_sct<true>,
                                           sdsl::sd_vector<>,
                                           RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQIlcpDSLPCore = dret::rmq::IlcpCore<TStorage,
                                            dret::Alphabet<>::int_width,
                                            sdsl::sd_vector<>,
                                            sdsl::rmq_succinct_sct<true>,
                                            sdsl::sd_vector<>,
                                            RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPCore = dret::rmq::CilcpCore<TStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::sd_vector<>,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQCilcpDSLPCore = dret::rmq::CilcpCore<TStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::sd_vector<>,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSCore = dret::rmq::SadaCore<TStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              RMQGetDocSLP_NS<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSCore = dret::rmq::IlcpCore<TStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::sd_vector<>,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              RMQGetDocSLP_NS<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSCore = dret::rmq::CilcpCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSRawCore = dret::rmq::SadaCore<TStorage,
                                                 dret::Alphabet<>::int_width,
                                                 sdsl::rmq_succinct_sct<true>,
                                                 sdsl::sd_vector<>,
                                                 RMQGetDocSLP_NS_Raw<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSRawCore = dret::rmq::IlcpCore<TStorage,
                                                 dret::Alphabet<>::int_width,
                                                 sdsl::sd_vector<>,
                                                 sdsl::rmq_succinct_sct<true>,
                                                 sdsl::sd_vector<>,
                                                 RMQGetDocSLP_NS_Raw<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSRawCore = dret::rmq::CilcpCore<TStorage,
                                                   dret::Alphabet<>::int_width,
                                                   sdsl::sd_vector<>,
                                                   sdsl::rmq_succinct_sct<true>,
                                                   sdsl::sd_vector<>,
                                                   RMQGetDocSLP_NS_Raw<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSDVCore = dret::rmq::SadaCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_DV<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSDVCore = dret::rmq::IlcpCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_DV<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSDVCore = dret::rmq::CilcpCore<TStorage,
                                                  dret::Alphabet<>::int_width,
                                                  sdsl::sd_vector<>,
                                                  sdsl::rmq_succinct_sct<true>,
                                                  sdsl::sd_vector<>,
                                                  RMQGetDocSLP_NS_DV<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSVVCore = dret::rmq::SadaCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_VV<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSVVCore = dret::rmq::IlcpCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_VV<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSVVCore = dret::rmq::CilcpCore<TStorage,
                                                  dret::Alphabet<>::int_width,
                                                  sdsl::sd_vector<>,
                                                  sdsl::rmq_succinct_sct<true>,
                                                  sdsl::sd_vector<>,
                                                  RMQGetDocSLP_NS_VV<TStorage>>;

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

template <typename TStorage>
using RMQSadaSLPNSIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQSadaSLPNSCore<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQIlcpSLPNSCore<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQCilcpSLPNSCore<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSRawIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQSadaSLPNSRawCore<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSRawIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQIlcpSLPNSRawCore<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSRawIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQCilcpSLPNSRawCore<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSDVIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQSadaSLPNSDVCore<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSDVIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQIlcpSLPNSDVCore<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSDVIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQCilcpSLPNSDVCore<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSVVIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQSadaSLPNSVVCore<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSVVIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQIlcpSLPNSVVCore<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSVVIndex =
    dret::rmq::DocListIdxRMQ<TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, RMQCilcpSLPNSVVCore<TStorage>>;

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
    RMQSadaSLPNSIndex<dret::GenericStorage>,
    RMQIlcpSLPNSIndex<dret::GenericStorage>,
    RMQCilcpSLPNSIndex<dret::GenericStorage>,
    dret::DocListIdxSLP<dret::GenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                        BareSLP_Raw>,
    RMQSadaSLPNSRawIndex<dret::GenericStorage>,
    RMQIlcpSLPNSRawIndex<dret::GenericStorage>,
    RMQCilcpSLPNSRawIndex<dret::GenericStorage>,
    dret::DocListIdxSLP<dret::GenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                        BareSLP_DV>,
    dret::DocListIdxSLP<dret::GenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                        BareSLP_VV>,
    RMQSadaSLPNSDVIndex<dret::GenericStorage>,
    RMQIlcpSLPNSDVIndex<dret::GenericStorage>,
    RMQCilcpSLPNSDVIndex<dret::GenericStorage>,
    RMQSadaSLPNSVVIndex<dret::GenericStorage>,
    RMQIlcpSLPNSVVIndex<dret::GenericStorage>,
    RMQCilcpSLPNSVVIndex<dret::GenericStorage>,
    RMQSadaDSLPIndex<dret::GenericStorage>,
    RMQIlcpDSLPIndex<dret::GenericStorage>,
    RMQCilcpDSLPIndex<dret::GenericStorage>,
    dret::pdl::DocListIdxPDLPlain<dret::GenericStorage>,
    dret::pdl::DocListIdxPDLRP<dret::GenericStorage>,
    dret::pdl::DocListIdxPDLBC<dret::GenericStorage>>;

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
    RMQSadaSLPNSIndex<ExternalGenericStorage>,
    RMQIlcpSLPNSIndex<ExternalGenericStorage>,
    RMQCilcpSLPNSIndex<ExternalGenericStorage>,
    dret::DocListIdxSLP<ExternalGenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                        BareSLP_Raw>,
    RMQSadaSLPNSRawIndex<ExternalGenericStorage>,
    RMQIlcpSLPNSRawIndex<ExternalGenericStorage>,
    RMQCilcpSLPNSRawIndex<ExternalGenericStorage>,
    dret::DocListIdxSLP<ExternalGenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                        BareSLP_DV>,
    dret::DocListIdxSLP<ExternalGenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                        BareSLP_VV>,
    RMQSadaSLPNSDVIndex<ExternalGenericStorage>,
    RMQIlcpSLPNSDVIndex<ExternalGenericStorage>,
    RMQCilcpSLPNSDVIndex<ExternalGenericStorage>,
    RMQSadaSLPNSVVIndex<ExternalGenericStorage>,
    RMQIlcpSLPNSVVIndex<ExternalGenericStorage>,
    RMQCilcpSLPNSVVIndex<ExternalGenericStorage>,
    RMQSadaDSLPIndex<ExternalGenericStorage>,
    RMQIlcpDSLPIndex<ExternalGenericStorage>,
    RMQCilcpDSLPIndex<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLRP<ExternalGenericStorage>,
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage>>;

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

// ListDocsRMQScheme / MarkedReported unit tests live in
// doc_list_rmq_scheme_test.cpp — they don't depend on Config / Factory and
// are split out to keep this integration suite focused.
