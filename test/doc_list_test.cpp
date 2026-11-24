//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <limits>
#include <random>
#include <set>
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
#include "dret/doc_list/doc_list_gcda_backend.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"

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
using RMQSadaSLPCore = dret::rmq::SadaLCore<TStorage,
                                           dret::Alphabet<>::int_width,
                                           sdsl::rmq_succinct_sct<true>,
                                           sdsl::sd_vector<>,
                                           RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQSadaDSLPCore = dret::rmq::SadaLCore<TStorage,
                                            dret::Alphabet<>::int_width,
                                            sdsl::rmq_succinct_sct<true>,
                                            sdsl::sd_vector<>,
                                            RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPCore = dret::rmq::IlcpLCore<TStorage,
                                           dret::Alphabet<>::int_width,
                                           sdsl::sd_vector<>,
                                           sdsl::rmq_succinct_sct<true>,
                                           sdsl::sd_vector<>,
                                           RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQIlcpDSLPCore = dret::rmq::IlcpLCore<TStorage,
                                            dret::Alphabet<>::int_width,
                                            sdsl::sd_vector<>,
                                            sdsl::rmq_succinct_sct<true>,
                                            sdsl::sd_vector<>,
                                            RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPCore = dret::rmq::CilcpLCore<TStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::sd_vector<>,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             RMQGetDocSLP<TStorage>>;

template <typename TStorage>
using RMQCilcpDSLPCore = dret::rmq::CilcpLCore<TStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::sd_vector<>,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             RMQGetDocDSLP<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSCore = dret::rmq::SadaLCore<TStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              RMQGetDocSLP_NS<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSCore = dret::rmq::IlcpLCore<TStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::sd_vector<>,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              RMQGetDocSLP_NS<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSCore = dret::rmq::CilcpLCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSRawCore = dret::rmq::SadaLCore<TStorage,
                                                 dret::Alphabet<>::int_width,
                                                 sdsl::rmq_succinct_sct<true>,
                                                 sdsl::sd_vector<>,
                                                 RMQGetDocSLP_NS_Raw<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSRawCore = dret::rmq::IlcpLCore<TStorage,
                                                 dret::Alphabet<>::int_width,
                                                 sdsl::sd_vector<>,
                                                 sdsl::rmq_succinct_sct<true>,
                                                 sdsl::sd_vector<>,
                                                 RMQGetDocSLP_NS_Raw<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSRawCore = dret::rmq::CilcpLCore<TStorage,
                                                   dret::Alphabet<>::int_width,
                                                   sdsl::sd_vector<>,
                                                   sdsl::rmq_succinct_sct<true>,
                                                   sdsl::sd_vector<>,
                                                   RMQGetDocSLP_NS_Raw<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSDVCore = dret::rmq::SadaLCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_DV<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSDVCore = dret::rmq::IlcpLCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_DV<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSDVCore = dret::rmq::CilcpLCore<TStorage,
                                                  dret::Alphabet<>::int_width,
                                                  sdsl::sd_vector<>,
                                                  sdsl::rmq_succinct_sct<true>,
                                                  sdsl::sd_vector<>,
                                                  RMQGetDocSLP_NS_DV<TStorage>>;

template <typename TStorage>
using RMQSadaSLPNSVVCore = dret::rmq::SadaLCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_VV<TStorage>>;

template <typename TStorage>
using RMQIlcpSLPNSVVCore = dret::rmq::IlcpLCore<TStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                RMQGetDocSLP_NS_VV<TStorage>>;

template <typename TStorage>
using RMQCilcpSLPNSVVCore = dret::rmq::CilcpLCore<TStorage,
                                                  dret::Alphabet<>::int_width,
                                                  sdsl::sd_vector<>,
                                                  sdsl::rmq_succinct_sct<true>,
                                                  sdsl::sd_vector<>,
                                                  RMQGetDocSLP_NS_VV<TStorage>>;

// RMQ × RLCSA cores are intentionally NOT defined here: RLCSA's GSA-style SA
// ordering doesn't match dret's SA ordering, so SADA/ILCP's RMQ recursion
// would mis-stop on the swapped doc IDs. PDL × RLCSA below works (set-based).
// See include/dret/rmq/rmq_get_doc_rlcsa.h.

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

// PDL codec aliases parameterised on the get-doc backend, so the type lists can
// exercise the SLP / bare-SLP / differential get-doc paths (not just the DA
// default). The bare path (PDLGetDocsSLP_NS) is the newest; SLP / DSLP were
// previously only covered by the benchmark.
template <typename TStorage, typename TGetDocs>
using PDLPlainBk = dret::pdl::DocListIdxPDLPlain<
    TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, TGetDocs>;
template <typename TStorage, typename TGetDocs>
using PDLRPBk = dret::pdl::DocListIdxPDLRP<
    TStorage, dret::Alphabet<>, RMQCountIdx<TStorage>, TGetDocs>;

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
                               grammar::CombinedSLPWithUnitCover<
                                   grammar::CombinedSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
                                                        grammar::SampledSLP<>,
                                                        sdsl::int_vector<>>>>,
    dret::DocListIdxSLP<>,
    dret::dgcda::DocListIdxDGCDA<>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<sdsl::int_vector<>>>>>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<sdsl::int_vector<>>>>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaLCore<dret::GenericStorage>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpLCore<dret::GenericStorage>>,
    dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpLCore<dret::GenericStorage>>,
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
    dret::DocListIdxSLP<dret::GenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                        dret::DifferentialSLP<>>,  // bare-diff (iv, compressed base)
    dret::DocListIdxSLP<dret::GenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                        dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>, sdsl::dac_vector<>,
                                              sdsl::dac_vector<>, sdsl::dac_vector<>>>,  // bare-diff (dv)
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
    dret::pdl::DocListIdxPDLBC<dret::GenericStorage>,
    PDLPlainBk<dret::GenericStorage, dret::pdl::PDLGetDocsSLP<dret::GenericStorage>>,
    PDLPlainBk<dret::GenericStorage, dret::pdl::PDLGetDocsSLP_NS<dret::GenericStorage>>,
    PDLPlainBk<dret::GenericStorage, dret::pdl::PDLGetDocsDSLP<dret::GenericStorage>>,
    PDLRPBk<dret::GenericStorage, dret::pdl::PDLGetDocsSLP<dret::GenericStorage>>,
    PDLRPBk<dret::GenericStorage, dret::pdl::PDLGetDocsSLP_NS<dret::GenericStorage>>,
    PDLRPBk<dret::GenericStorage, dret::pdl::PDLGetDocsDSLP<dret::GenericStorage>>,
    // RMQ × RLCSA omitted (search-type incompatibility — see DocListIndexSearchTypes).
    // PDL × RLCSA still exercises the RLCSA construct+load path.
    PDLPlainBk<dret::GenericStorage, dret::pdl::PDLGetDocsRLCSA<dret::GenericStorage>>>;

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

// GCDA with plain document sets -- the sorted lists bit-packed, no grammar over
// them (GCDASetsCodec::Plain in the benchmark) -- over the default sampled-cached
// grammar and over the differential one.
using PlainSets = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;
template <typename S>
using GCDAPlainSets = dret::gcda::DocListIdxGCDA<
    S, dret::Alphabet<>, sri::RIndexCount<S, dret::Alphabet<>>,
    grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>, grammar::SampledSLP<>,
                      grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
    PlainSets>;
template <typename S>
using DGCDAPlainSets = dret::dgcda::DocListIdxDGCDA<
    S, dret::Alphabet<>, sri::RIndexCount<S, dret::Alphabet<>>, dret::DifferentialLightSLP<>, PlainSets>;

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
    GCDAPlainSets<ExternalGenericStorage>,
    DGCDAPlainSets<ExternalGenericStorage>,
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
                               grammar::CombinedSLPWithUnitCover<
                                   grammar::CombinedSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
                                                        grammar::SampledSLP<>,
                                                        sdsl::int_vector<>>>>,
    dret::DocListIdxSLP<ExternalGenericStorage>,
    dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage>,
    dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<sdsl::int_vector<>>>>>,
    dret::dgcda::DocListIdxDGCDA<dret::GenericStorage,
                                 dret::Alphabet<>,
                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                 dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<sdsl::int_vector<>>>>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaLCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpLCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                             dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpLCore<ExternalGenericStorage>>,
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
    dret::DocListIdxSLP<ExternalGenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                        dret::DifferentialSLP<>>,  // bare-diff (iv, compressed base)
    dret::DocListIdxSLP<ExternalGenericStorage,
                        dret::Alphabet<>,
                        sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                        dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>, sdsl::dac_vector<>,
                                              sdsl::dac_vector<>, sdsl::dac_vector<>>>,  // bare-diff (dv)
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
    dret::pdl::DocListIdxPDLBC<ExternalGenericStorage>,
    PDLPlainBk<ExternalGenericStorage, dret::pdl::PDLGetDocsSLP<ExternalGenericStorage>>,
    PDLPlainBk<ExternalGenericStorage, dret::pdl::PDLGetDocsSLP_NS<ExternalGenericStorage>>,
    PDLPlainBk<ExternalGenericStorage, dret::pdl::PDLGetDocsDSLP<ExternalGenericStorage>>,
    PDLRPBk<ExternalGenericStorage, dret::pdl::PDLGetDocsSLP<ExternalGenericStorage>>,
    PDLRPBk<ExternalGenericStorage, dret::pdl::PDLGetDocsSLP_NS<ExternalGenericStorage>>,
    PDLRPBk<ExternalGenericStorage, dret::pdl::PDLGetDocsDSLP<ExternalGenericStorage>>,
    // RMQ × RLCSA is NOT in the search-type list: RLCSA's GSA-style SA
    // ordering doesn't match dret's, and SADA/ILCP miss docs on real data.
    // (Small test data doesn't trip the bug because patterns avoid the
    // "ambiguous zones" where multiple doc-end suffixes share a prefix.)
    // The factory falls RMQ × RLCSA through to DA; see rmq_get_doc_rlcsa.h.
    // PDL × RLCSA is set-based and works correctly.
    PDLPlainBk<ExternalGenericStorage, dret::pdl::PDLGetDocsRLCSA<ExternalGenericStorage>>,
    PDLRPBk<ExternalGenericStorage, dret::pdl::PDLGetDocsRLCSA<ExternalGenericStorage>>>;

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

// Modification times of a cached object's component files; every one must exist.
std::vector<std::filesystem::file_time_type> ComponentTimes(const std::vector<dret::Component>& t_components,
                                                            const sdsl::cache_config& t_config) {
  std::vector<std::filesystem::file_time_type> times;
  for (const auto& c : t_components) {
    const auto path = dret::ComponentFile(c, t_config);
    EXPECT_TRUE(std::filesystem::exists(path)) << path;
    times.push_back(std::filesystem::exists(path) ? std::filesystem::last_write_time(path)
                                                  : std::filesystem::file_time_type{});
  }
  return times;
}

// Number of cache files whose name starts with t_prefix.
std::size_t CountCacheFiles(const sdsl::cache_config& t_config, const std::string& t_prefix) {
  std::size_t n = 0;
  for (const auto& e : std::filesystem::directory_iterator(t_config.dir))
    n += e.path().filename().string().starts_with(t_prefix);
  return n;
}

// Number of files of the component t_key, one per encoding: "<key>_<type hash>_...".
std::size_t CountComponentFiles(const sdsl::cache_config& t_config, const std::string& t_key) {
  std::size_t n = 0;
  for (const auto& e : std::filesystem::directory_iterator(t_config.dir)) {
    const auto name = e.path().filename().string();
    n += name.starts_with(t_key + "_") && name.size() > t_key.size() + 1 && std::isdigit(name[t_key.size() + 1]);
  }
  return n;
}

const std::vector<std::pair<std::string, std::vector<std::size_t>>> kReuseExpected = {
    {"TAT", {0}}, {"LAT", {1}}, {"LAL", {2}}, {"TA", {0, 1}}, {"LA", {1, 2}}, {"A", {0, 1, 2}}, {"TAL", {}},
};

template <typename TIndex>
void ExpectReuseResults(const TIndex& t_index, const std::string& t_label) {
  for (const auto& [pattern, docs] : kReuseExpected) {
    DocListResultVector result;
    t_index.Search(pattern, std::ref(result));
    result();
    EXPECT_THAT(result, testing::ElementsAreArray(docs)) << t_label << " " << pattern;
  }
}

TEST_F(RMQSLPCacheReuseTest, rmq_slp_reuses_gcda_slp_cache) {
  using GetDocSLP = RMQGetDocSLP<ExternalGenericStorage>;
  using TSLP = typename GetDocSLP::SLP;
  using RMQSadaSLP = RMQSadaSLPIndex<ExternalGenericStorage>;

  dret::gcda::DocListIdxGCDA<ExternalGenericStorage> gcda(std::ref(storage_), 512, 4);
  construct(gcda, config_);

  const auto components = dret::CacheComponents(TSLP{}, config_.keys, dret::SampledTreeCell{512, 4});
  const auto before = ComponentTimes(components, config_);

  RMQSadaSLP rmq_sada_slp(std::ref(storage_));
  construct(rmq_sada_slp, config_);

  EXPECT_EQ(ComponentTimes(components, config_), before);
}

// The RMQ backend for GCDA's differential grammar without the sampled tree: a
// bare DifferentialSLP whose block size is the sample spacing. Each block size
// answers correctly and stores its own samples; the grammar, roots and span
// sums are stored once for all of them.
TEST_F(RMQSLPCacheReuseTest, rmq_bare_diff_sweeps_block_size) {
  using S = ExternalGenericStorage;
  using TDiff = dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using TGetDoc = dret::rmq::GetDocSLP_NS<S, dret::Alphabet<>::int_width, TDiff>;
  using TCount = RMQCountIdx<S>;
  using SadaL = dret::rmq::SadaLCore<S, dret::Alphabet<>::int_width, sdsl::rmq_succinct_sct<true>,
                                     sdsl::sd_vector<>, TGetDoc>;
  using IlcpL = dret::rmq::IlcpLCore<S, dret::Alphabet<>::int_width, sdsl::sd_vector<>,
                                     sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
  using CilcpL = dret::rmq::CilcpLCore<S, dret::Alphabet<>::int_width, sdsl::sd_vector<>,
                                       sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;

  auto check = [&]<typename TCore>(std::uint32_t bs, const char* core_name) {
    TCore core(std::ref(storage_), bs, 0.0f);
    dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, TCount, TCore> index(std::ref(storage_), core);
    construct(index, config_);
    index.load(config_);
    ExpectReuseResults(index, std::string(core_name) + " bs=" + std::to_string(bs));
  };

  for (std::uint32_t bs : {2u, 4u, 512u}) {
    check.template operator()<SadaL>(bs, "SADA-L");
    check.template operator()<IlcpL>(bs, "ILCP-L");
    check.template operator()<CilcpL>(bs, "CILCP-L");
    EXPECT_TRUE(dret::ComponentsExist(dret::CacheComponents(TDiff{}, config_.keys, bs), config_)) << bs;
  }
  EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, dret::conf::kDaDiffGrammar)), 1u);
  EXPECT_EQ(CountCacheFiles(config_, "spc"), 3u);
}

// At the GCDA-nolists spacing the RMQ backend loads the components the
// GCDA-nolists differential index (DocListIdxSLP) built, not rebuild its own.
TEST_F(RMQSLPCacheReuseTest, rmq_bare_diff_reuses_gcda_nolists_grammar) {
  using S = ExternalGenericStorage;
  using TDiff = dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using TGetDoc = dret::rmq::GetDocSLP_NS<S, dret::Alphabet<>::int_width, TDiff>;
  using CilcpL = dret::rmq::CilcpLCore<S, dret::Alphabet<>::int_width, sdsl::sd_vector<>,
                                       sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;

  dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, TDiff> nolists(std::ref(storage_));
  construct(nolists, config_);

  const auto components = dret::CacheComponents(TDiff{}, config_.keys, dret::kDiffBlockSize);
  const auto before = ComponentTimes(components, config_);

  CilcpL core(std::ref(storage_), dret::kDiffBlockSize, 0.0f);
  dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, RMQCountIdx<S>, CilcpL> rmq(std::ref(storage_), core);
  construct(rmq, config_);

  EXPECT_EQ(ComponentTimes(components, config_), before);
  EXPECT_EQ(CountCacheFiles(config_, "spc"), 1u);
}

// A PDL-RP core cached with the pre-2026-10 codec (RPCodecWithLengths, whose
// lists carry a per-rule length array they never read) and the current one
// share the tree and the selection; the current lists answer the same and are
// no larger.
TEST_F(RMQSLPCacheReuseTest, pdl_rp_lists_without_lengths) {
  using S = ExternalGenericStorage;
  using Legacy = dret::pdl::DocListIdxPDL<S, dret::Alphabet<>, sri::RIndexCount<S, dret::Alphabet<>>,
                                          dret::pdl::PDLGetDocsDA<S>, dret::pdl::RPCodecWithLengths>;
  using Current = dret::pdl::DocListIdxPDLRP<S>;

  Legacy legacy(std::ref(storage_), 2, 2.0f);
  construct(legacy, config_);

  Current current(std::ref(storage_), 2, 2.0f);
  construct(current, config_);

  EXPECT_LE(sdsl::size_in_bytes(current), sdsl::size_in_bytes(legacy));
  ExpectReuseResults(current, "PDL-RP");
  EXPECT_EQ(CountCacheFiles(config_, "blk2_" + dret::KeyName(config_.keys, dret::conf::kPdlTree)), 1u);
  EXPECT_EQ(CountCacheFiles(config_, "blk2-sf2-occw_" + dret::KeyName(config_.keys, dret::conf::kPdlSelection)), 1u);
}

// GCDA-nolists' differential index sweeps the same sample spacing as the RMQ
// backend: each spacing answers correctly and stores its own samples, and all
// of them share one grammar, one roots and one span-sums component.
TEST_F(RMQSLPCacheReuseTest, gcda_nolists_diff_sweeps_spacing) {
  using S = ExternalGenericStorage;
  using TDiff = dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using Index = dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, TDiff>;
  std::vector<std::size_t> sizes;
  for (std::uint32_t spacing : {2u, 4u, 512u}) {
    Index index(std::ref(storage_), spacing);
    construct(index, config_);
    sizes.push_back(sdsl::size_in_bytes(index));
    ExpectReuseResults(index, "spacing " + std::to_string(spacing));
    EXPECT_TRUE(dret::ComponentsExist(dret::CacheComponents(TDiff{}, config_.keys, spacing), config_)) << spacing;
  }
  for (auto name : {dret::conf::kDaDiffGrammar, dret::conf::kDaDiffRoots, dret::conf::kDaDiffSpanSums})
    EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, name)), 1u) << name;
  EXPECT_EQ(CountCacheFiles(config_, "spc"), 3u);
  // The spacing is real: a denser sampling stores more samples.
  EXPECT_GT(sizes[0], sizes[2]);
}

// The differential container variants share the grammar; each container keeps
// its own encoding of the roots, span sums and samples, and the int_vector
// roots (the sequence a new spacing restarts from) are stored once.
TEST_F(RMQSLPCacheReuseTest, gcda_nolists_diff_containers_share_the_grammar) {
  using S = ExternalGenericStorage;
  using TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
  using TDiffIV = dret::DifferentialSLP<TSLP>;
  using TDiffDV = dret::DifferentialSLP<TSLP, sdsl::dac_vector<>, sdsl::dac_vector<>, sdsl::dac_vector<>>;
  dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, TDiffDV> dv(std::ref(storage_), 4);
  construct(dv, config_);
  dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, TDiffIV> iv(std::ref(storage_), 4);
  construct(iv, config_);
  ExpectReuseResults(dv, "diff-dv");
  ExpectReuseResults(iv, "diff-iv");
  EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, dret::conf::kDaDiffGrammar)), 1u);
  EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, dret::conf::kDaDiffRoots)), 2u);
}

// Every GCDA representation of a cell, and GCDA-differential, read one sampled
// tree; the on-demand representation and GCDA-nolists (plain) one CNF grammar.
TEST_F(RMQSLPCacheReuseTest, gcda_representations_share_components) {
  using S = ExternalGenericStorage;
  using SLPCNF = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
  using Combined = grammar::CombinedSLPWithUnitCover<grammar::CombinedSLP<SLPCNF, grammar::SampledSLP<>, sdsl::int_vector<>>>;
  using Sets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>, true,
                                 grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using SetsPlain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;
  using RC = sri::RIndexCount<S, dret::Alphabet<>>;
  const dret::SampledTreeCell cell{2, 2};

  dret::gcda::DocListIdxGCDA<S> light(std::ref(storage_), cell.block_size, cell.storing_factor);
  EXPECT_NO_THROW(construct(light, config_)) << "light";
  dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, Combined, Sets> combined(std::ref(storage_), 2, 2);
  EXPECT_NO_THROW(construct(combined, config_)) << "combined";
  dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, grammar::CompactBPSLP<>, SetsPlain> bp(std::ref(storage_), 2, 2);
  EXPECT_NO_THROW(construct(bp, config_)) << "bp";
  dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, grammar::CompactLOUDSSLP<>, Sets> louds(std::ref(storage_), 2, 2);
  EXPECT_NO_THROW(construct(louds, config_)) << "louds";
  dret::dgcda::DocListIdxDGCDA<S> dgcda(std::ref(storage_), 2, 2);
  EXPECT_NO_THROW(construct(dgcda, config_)) << "dgcda";
  dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, SLPCNF> nolists(std::ref(storage_));
  EXPECT_NO_THROW(construct(nolists, config_)) << "nolists";

  ExpectReuseResults(light, "cached");
  ExpectReuseResults(combined, "ondemand");
  ExpectReuseResults(bp, "bp");
  ExpectReuseResults(louds, "louds");
  ExpectReuseResults(dgcda, "differential");
  ExpectReuseResults(nolists, "nolists");

  const auto tree = dret::CellKey(config_.keys, dret::conf::kDaSampledTree, cell);
  EXPECT_EQ(CountComponentFiles(config_, tree), 1u);
  // The CNF grammar: its raw encoding (built from irepair's output) and the
  // bit-compressed one, which on-demand GCDA and GCDA-nolists share.
  EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, dret::conf::kDaCnfGrammar)), 2u);
  // The BP and LOUDS encodings of it, one each for the whole collection.
  EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, dret::conf::kDaCnfGrammarBP)), 1u);
  EXPECT_EQ(CountComponentFiles(config_, dret::KeyName(config_.keys, dret::conf::kDaCnfGrammarLOUDS)), 1u);
  EXPECT_TRUE(std::filesystem::exists(dret::ComponentFile(dret::CacheComponents(SLPCNF{}, config_.keys).front(), config_)));
  // One list per codec and encoding: raw and bit-packed plain lists, and the
  // raw and final Re-Pair lists. GCDA and GCDA-differential share all of them.
  EXPECT_EQ(CountComponentFiles(config_, dret::CellKey(config_.keys, dret::conf::kDaNodeDocListsPlain, cell)), 2u);
  EXPECT_EQ(CountComponentFiles(config_, dret::CellKey(config_.keys, dret::conf::kDaNodeDocListsRP, cell)), 2u);
}

// A leaves-only PDL core does not read the storing factor: the cores of every
// storing factor are one set of files. An occurrence-weighted core gets its
// own selection and lists per storing factor, and all of them share the tree.
TEST_F(RMQSLPCacheReuseTest, pdl_components_follow_the_parameters_read) {
  using S = ExternalGenericStorage;
  using PDL = dret::pdl::DocListIdxPDLPlain<S>;
  for (float sf : {2.0f, 4.0f}) {
    for (auto policy : {dret::pdl::StoragePolicy::LeavesOnly, dret::pdl::StoragePolicy::OccurrenceWeighted}) {
      PDL pdl(std::ref(storage_), 2, sf, policy);
      construct(pdl, config_);
      ExpectReuseResults(pdl, "sf " + std::to_string(sf));
    }
  }
  EXPECT_EQ(CountCacheFiles(config_, "blk2_" + dret::KeyName(config_.keys, dret::conf::kPdlTree)), 1u);
  EXPECT_EQ(CountCacheFiles(config_, "blk2-leaves_"), 2u);  // selection, lists
  EXPECT_EQ(CountCacheFiles(config_, "blk2-sf2-occw_"), 2u);
  EXPECT_EQ(CountCacheFiles(config_, "blk2-sf4-occw_"), 2u);
}

TEST(PolicyReadsStoringFactorTest, only_occurrence_weighted) {
  EXPECT_TRUE(dret::pdl::PolicyReadsStoringFactor(dret::pdl::StoragePolicy::OccurrenceWeighted));
  EXPECT_FALSE(dret::pdl::PolicyReadsStoringFactor(dret::pdl::StoragePolicy::LeavesOnly));
  EXPECT_FALSE(dret::pdl::PolicyReadsStoringFactor(dret::pdl::StoragePolicy::StoreAllInternal));
}

// The file names: each prefix carries the parameters it names.
TEST(PrefixedKeyTest, names_carry_their_parameters) {
  const auto keys = dret::Keys<8>().keys;
  EXPECT_EQ(dret::CellKey(keys, dret::conf::kDaSampledTree, {512, 4}), "blk512-sf4_da_sampled_tree");
  EXPECT_EQ(dret::CellKey(keys, dret::conf::kDaNodeDocListsRP, {1024, 32}), "blk1024-sf32_da_node_doclists_rp");
  EXPECT_EQ(dret::PrefixedKey(keys, dret::conf::kSpc, "da_diff_samples", 128u), "spc128_da_diff_samples");
  // beta = infinity names itself, for GCDA's cells and PDL's selections.
  EXPECT_EQ(dret::CellKey(keys, dret::conf::kDaSampledTree, {512, dret::kInfiniteStoringFactor}),
            "blk512-sfinf_da_sampled_tree");
  EXPECT_EQ(dret::pdl::PdlSelectionKey(keys, dret::conf::kPdlSelection, 256, dret::kInfiniteStoringFactor,
                                       dret::pdl::StoragePolicy::OccurrenceWeighted),
            "blk256-sfinf-occw_pdl_selection");
  EXPECT_EQ(dret::PrefixedKey(keys, dret::conf::kBlk, "pdl_tree", 256u), "blk256_pdl_tree");
  using dret::pdl::StoragePolicy;
  EXPECT_EQ(dret::pdl::PdlSelectionKey(keys, dret::conf::kPdlSelection, 256, 8, StoragePolicy::OccurrenceWeighted),
            "blk256-sf8-occw_pdl_selection");
  EXPECT_EQ(dret::pdl::PdlSelectionKey(keys, dret::conf::kPdlDocListsRP, 256, 8, StoragePolicy::LeavesOnly),
            "blk256-leaves_pdl_doclists_rp");
}

// A component that exists with another length means its key misses a parameter
// of its construction: storing over it throws instead of overwriting. One of the
// same length is kept as it is (see StoreComponents).
TEST_F(RMQSLPCacheReuseTest, store_components_rejects_a_different_component) {
  sdsl::int_vector<> v{1, 2, 3};
  const std::vector<dret::Component> one{{"probe", dret::TypeHash<sdsl::int_vector<>>(), dret::SerializedSize(v)}};
  dret::StoreComponents(v, one, config_);
  EXPECT_NO_THROW(dret::StoreComponents(v, one, config_));
  sdsl::int_vector<> w{1, 2, 3, 4};
  EXPECT_THROW(dret::StoreComponents(w, {{"probe", one[0].type, dret::SerializedSize(w)}}, config_),
               std::logic_error);
  sdsl::int_vector<> loaded;
  ASSERT_TRUE(dret::LoadComponents(loaded, one, config_));
  EXPECT_EQ(loaded, v);
  // A layout that does not cover the whole serialization is rejected too.
  EXPECT_THROW(dret::StoreComponents(v, {{"probe2", one[0].type, one[0].bytes - 1}}, config_), std::logic_error);
}

TEST_F(RMQSLPCacheReuseTest, rmq_dslp_reuses_dgcda_dslp_cache) {
  using GetDocDSLP = RMQGetDocDSLP<ExternalGenericStorage>;
  using TDSLP = typename GetDocDSLP::DSLP;
  using RMQSadaDSLP = RMQSadaDSLPIndex<ExternalGenericStorage>;

  dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage> dgcda(std::ref(storage_), 512, 4);
  construct(dgcda, config_);

  const auto components = dret::CacheComponents(TDSLP{}, config_.keys, dret::SampledTreeCell{512, 4});
  const auto before = ComponentTimes(components, config_);

  RMQSadaDSLP rmq_sada_dslp(std::ref(storage_));
  construct(rmq_sada_dslp, config_);

  EXPECT_EQ(ComponentTimes(components, config_), before);
}

TEST_F(RMQSLPCacheReuseTest, pdl_rlcsa_sidecar_not_rebuilt) {
  // The RLCSA sidecar is constructed by GetDocRLCSA via PDLRawRangePolicy.
  // Test via PDL × RLCSA — the only supported consumer (RMQ × RLCSA is
  // disabled at the factory level due to SA-ordering incompatibility; see
  // rmq_get_doc_rlcsa.h header comment).
  using PDLRLCSA = PDLPlainBk<ExternalGenericStorage, dret::pdl::PDLGetDocsRLCSA<ExternalGenericStorage>>;

  PDLRLCSA pdl1(std::ref(storage_));
  construct(pdl1, config_);

  const auto base = sdsl::cache_file_name(
      config_.keys[dret::conf::kRLCSA].get<std::string>(), config_);
  const auto array_path = base + ".rlcsa.array";
  ASSERT_TRUE(std::filesystem::exists(array_path));
  const auto before_time = std::filesystem::last_write_time(array_path);
  const auto before_size = std::filesystem::file_size(array_path);

  PDLRLCSA pdl2(std::ref(storage_));
  construct(pdl2, config_);

  EXPECT_EQ(std::filesystem::file_size(array_path), before_size);
  EXPECT_EQ(std::filesystem::last_write_time(array_path), before_time);
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


// Every family built over a grid of its parameters in ONE cache directory, on
// a repetitive collection large enough for the cells' structures to differ.
// Shared components are keyed by the parameters their construction reads. Each
// cell is also built alone, in a private cache directory: the two must be the
// same index (same size, same answers). A key missing a parameter makes a cell
// read another cell's component, and that changes the size even where it cannot
// change an answer (a PDL selection, say).
class ComponentSweepTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    std::mt19937 rng(20261004);
    std::uniform_int_distribution<int> base(0, 3), pos(0, 239), coin(0, 9);
    std::string ref;
    for (int i = 0; i < 240; ++i) ref += "ACGT"[base(rng)];
    for (int d = 0; d < 24; ++d) {
      auto doc = ref;
      for (int k = 0; k < 6; ++k) doc[pos(rng)] = "ACGT"[base(rng)];
      if (coin(rng) < 3) doc = doc.substr(0, 120 + pos(rng) / 2);
      docs_.push_back(doc);
    }
    std::string data;
    for (const auto& d : docs_) { data += d; data += '\1'; }
    Init(data);
    pristine_ = config_;

    std::set<std::string> pats;
    for (const auto& doc : docs_)
      for (std::size_t i = 0; i + 8 <= doc.size(); i += 7)
        for (std::size_t m : {1u, 3u, 6u, 8u}) pats.insert(doc.substr(i, m));
    patterns_.assign(pats.begin(), pats.end());
  }

  template <typename TIndex>
  void CheckAnswers(const TIndex& t_index, const std::string& t_label) {
    std::size_t wrong = 0;
    for (const auto& p : patterns_) {
      std::vector<std::size_t> expected;
      for (std::size_t d = 0; d < docs_.size(); ++d)
        if (docs_[d].find(p) != std::string::npos) expected.push_back(d);
      DocListResultVector result;
      t_index.Search(p, std::ref(result));
      result();
      wrong += !std::equal(result.begin(), result.end(), expected.begin(), expected.end());
    }
    EXPECT_EQ(wrong, 0u) << t_label;
  }

  // t_build(config, storage) constructs one index and returns it. Built in the
  // shared directory and in a private one; both must answer right and match.
  template <typename TBuild>
  void Cell(const std::string& t_label, TBuild&& t_build) {
    sri::GenericStorage shared_storage;
    std::size_t shared_size = 0;
    EXPECT_NO_THROW({
      auto shared = t_build(config_, shared_storage);
      CheckAnswers(shared, t_label);
      shared_size = sdsl::size_in_bytes(shared);
    }) << t_label;

    auto fresh_config = pristine_;
    fresh_config.dir = (tmp_dir_ / ("fresh" + std::to_string(n_fresh_++))).string();
    std::filesystem::create_directories(fresh_config.dir);
    sri::GenericStorage fresh_storage;
    auto fresh = t_build(fresh_config, fresh_storage);
    EXPECT_EQ(shared_size, sdsl::size_in_bytes(fresh)) << t_label << ": the shared cache gave another index";
  }

  dret::Config pristine_;
  std::size_t n_fresh_ = 0;
  std::vector<std::string> docs_;
  std::vector<std::string> patterns_;
};

TEST_F(ComponentSweepTest, gcda_family_grid) {
  using S = ExternalGenericStorage;
  using RC = sri::RIndexCount<S, dret::Alphabet<>>;
  using SLPCNF = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
  using Light = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>, grammar::SampledSLP<>,
                                  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using Combined = grammar::CombinedSLPWithUnitCover<grammar::CombinedSLP<SLPCNF, grammar::SampledSLP<>, sdsl::int_vector<>>>;
  using SetsRP = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>, true,
                                   grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using SetsPlain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;
  for (std::uint32_t b : {4u, 16u, 64u}) {
    for (float sf : {1.0f, 4.0f, dret::kInfiniteStoringFactor}) {
      const auto label = " b=" + std::to_string(b) + " sf=" + std::to_string(sf);
      auto run = [&]<typename TIndex>(const std::string& name) {
        Cell(name + label, [&](dret::Config& cfg, sri::GenericStorage& st) {
          TIndex index(std::ref(st), b, sf);
          construct(index, cfg);
          return index;
        });
      };
      run.template operator()<dret::gcda::DocListIdxGCDA<S>>("cached");
      run.template operator()<dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, Light, SetsPlain>>("cached-plain");
      run.template operator()<dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, Combined, SetsRP>>("ondemand");
      run.template operator()<dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, grammar::CompactBPSLP<>, SetsRP>>("bp");
      run.template operator()<dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, grammar::CompactLOUDSSLP<>, SetsPlain>>("louds");
      run.template operator()<dret::dgcda::DocListIdxDGCDA<S>>("differential");
    }
  }
}

TEST_F(ComponentSweepTest, nolists_and_rmq_backends_grid) {
  using S = ExternalGenericStorage;
  using TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
  using TDiff = dret::DifferentialSLP<TSLP>;
  using TDiffDV = dret::DifferentialSLP<TSLP, sdsl::dac_vector<>, sdsl::dac_vector<>, sdsl::dac_vector<>>;
  auto nolists = [&]<typename T>(const std::string& name, std::uint32_t spacing) {
    Cell(name + " spc=" + std::to_string(spacing), [&](dret::Config& cfg, sri::GenericStorage& st) {
      dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, T> index(std::ref(st), spacing);
      construct(index, cfg);
      return index;
    });
  };
  for (std::uint32_t spacing : {2u, 8u, 64u}) {
    nolists.template operator()<TDiff>("nolists-diff", spacing);
    nolists.template operator()<TDiffDV>("nolists-diff-dv", spacing);
  }
  nolists.template operator()<TSLP>("nolists-plain", dret::kDiffBlockSize);

  using TGetDocDiff = dret::rmq::GetDocSLP_NS<S, dret::Alphabet<>::int_width, TDiff>;
  using CilcpLDiff = dret::rmq::CilcpLCore<S, dret::Alphabet<>::int_width, sdsl::sd_vector<>,
                                           sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDocDiff>;
  using IlcpLSLP = dret::rmq::IlcpLCore<S, dret::Alphabet<>::int_width, sdsl::sd_vector<>,
                                        sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, RMQGetDocSLP<S>>;
  auto rmq = [&]<typename TCore>(const std::string& name, std::uint32_t b, float sf) {
    Cell(name + " b=" + std::to_string(b) + " sf=" + std::to_string(sf),
         [&](dret::Config& cfg, sri::GenericStorage& st) {
           TCore core(std::ref(st), b, sf);
           dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, RMQCountIdx<S>, TCore> index(std::ref(st), core);
           construct(index, cfg);
           index.load(cfg);
           return index;
         });
  };
  for (std::uint32_t b : {4u, 32u}) {
    rmq.template operator()<CilcpLDiff>("CILCP-L diff", b, 0.0f);
    for (float sf : {2.0f, 32.0f}) rmq.template operator()<IlcpLSLP>("ILCP-L cached", b, sf);
  }
}

TEST_F(ComponentSweepTest, pdl_grid) {
  using S = ExternalGenericStorage;
  using dret::pdl::StoragePolicy;
  for (std::uint32_t b : {4u, 32u}) {
    for (float sf : {1.0f, 8.0f, dret::kInfiniteStoringFactor}) {
      for (auto policy : {StoragePolicy::OccurrenceWeighted, StoragePolicy::LeavesOnly, StoragePolicy::ListWeighted}) {
        const auto label = " b=" + std::to_string(b) + " sf=" + std::to_string(sf) + " policy=" +
                           std::to_string(static_cast<int>(policy));
        auto run = [&]<typename TIndex>(const std::string& name) {
          Cell(name + label, [&](dret::Config& cfg, sri::GenericStorage& st) {
            TIndex index(std::ref(st), b, sf, policy);
            construct(index, cfg);
            return index;
          });
        };
        run.template operator()<dret::pdl::DocListIdxPDLPlain<S>>("PDL-plain");
        run.template operator()<dret::pdl::DocListIdxPDLRP<S>>("PDL-rp");
      }
    }
  }
}

// GCDA over a backend: the sampled tree and node lists of GCDA, the ends of a
// range expanded by the plain DA or SA-Phi. Same answers as brute force, and the
// shared cache gives the same index as a private one.
TEST_F(ComponentSweepTest, gcda_backend_grid) {
  using S = ExternalGenericStorage;
  using RC = sri::RIndexCount<S, dret::Alphabet<>>;
  using SetsPlain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;
  using SetsRP = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>, true,
                                   grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using DA = dret::pdl::PDLGetDocsDA<S, dret::Alphabet<>::int_width>;
  using SAPhi = dret::pdl::PDLGetDocsSAPhi_R<S, dret::Alphabet<>::int_width>;
  for (std::uint32_t b : {4u, 64u}) {
    for (float sf : {1.0f, 4.0f, dret::kInfiniteStoringFactor}) {
      const auto label = " b=" + std::to_string(b) + " sf=" + std::to_string(sf);
      auto run = [&]<typename TIndex>(const std::string& name) {
        Cell(name + label, [&](dret::Config& cfg, sri::GenericStorage& st) {
          TIndex index(std::ref(st), b, sf);
          construct(index, cfg);
          return index;
        });
      };
      run.template operator()<dret::gcda::DocListIdxGCDABackend<S, dret::Alphabet<>, RC, DA, SetsRP>>("GCDA over DA, rp");
      run.template operator()<dret::gcda::DocListIdxGCDABackend<S, dret::Alphabet<>, RC, SAPhi, SetsPlain>>("GCDA over SA-Phi, plain");
      run.template operator()<dret::gcda::DocListIdxGCDABackend<S, dret::Alphabet<>, RC, SAPhi, SetsRP>>("GCDA over SA-Phi, rp");
    }
  }
}
