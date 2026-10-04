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

// The RMQ backend for GCDA's differential grammar without the sampled tree: a
// bare DifferentialSLP whose block size is the sample spacing. Each block size
// must answer correctly and get its own cache entry, so the cells of a
// block-size sweep never load one another's grammar -- except at the GCDA-nolists
// block size, which keeps that index's plain kDSLPNS key (see the next test).
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

  const std::vector<std::pair<std::string, std::vector<std::size_t>>> expected = {
      {"TAT", {0}}, {"LAT", {1}}, {"LAL", {2}}, {"TA", {0, 1}},
      {"LA", {1, 2}}, {"A", {0, 1, 2}}, {"TAL", {}},
  };

  auto check = [&]<typename TCore>(std::uint32_t bs, const char* core_name) {
    TCore core(std::ref(storage_), bs, 0.0f);
    dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, TCount, TCore> index(std::ref(storage_), core);
    construct(index, config_);
    index.load(config_);
    for (const auto& [pattern, docs] : expected) {
      DocListResultVector result;
      index.Search(pattern, std::ref(result));
      result();
      EXPECT_THAT(result, testing::ElementsAreArray(docs)) << core_name << " bs=" << bs << " " << pattern;
    }
  };

  for (std::uint32_t bs : {2u, 4u, 512u}) {
    check.template operator()<SadaL>(bs, "SADA-L");
    check.template operator()<IlcpL>(bs, "ILCP-L");
    check.template operator()<CilcpL>(bs, "CILCP-L");

    const auto key = dret::DiffNoTreeCacheKey(config_.keys, bs);
    EXPECT_EQ(key.starts_with("bs"), bs != dret::kDiffBlockSize) << key;
    EXPECT_TRUE(std::filesystem::exists(sdsl::cache_file_name<TDiff>(key, config_))) << key;
  }

  // The RePair base grammar does not depend on the block size: one cached copy,
  // under the GCDA-nolists key, serves every block size.
  std::size_t n_grammars = 0;
  const auto dir = std::filesystem::path(
      sdsl::cache_file_name<TDiff>(dret::DiffNoTreeCacheKey(config_.keys, 4u), config_)).parent_path();
  for (const auto& e : std::filesystem::directory_iterator(dir))
    n_grammars += e.path().filename().string().find("_grammar_") != std::string::npos;
  EXPECT_EQ(n_grammars, 1u);
}

// At the GCDA-nolists block size the RMQ backend must load the grammar file that
// the GCDA-nolists differential index (DocListIdxSLP) built, not rebuild its own.
TEST_F(RMQSLPCacheReuseTest, rmq_bare_diff_reuses_gcda_nolists_grammar) {
  using S = ExternalGenericStorage;
  using TDiff = dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using TGetDoc = dret::rmq::GetDocSLP_NS<S, dret::Alphabet<>::int_width, TDiff>;
  using CilcpL = dret::rmq::CilcpLCore<S, dret::Alphabet<>::int_width, sdsl::sd_vector<>,
                                       sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;

  dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, TDiff> nolists(std::ref(storage_));
  construct(nolists, config_);

  const auto key = config_.keys[dret::conf::kDSLPNS].get<std::string>();
  const auto path = sdsl::cache_file_name<TDiff>(key, config_);
  ASSERT_TRUE(std::filesystem::exists(path));
  const auto before_time = std::filesystem::last_write_time(path);

  CilcpL core(std::ref(storage_), dret::kDiffBlockSize, 0.0f);
  dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, RMQCountIdx<S>, CilcpL> rmq(std::ref(storage_), core);
  construct(rmq, config_);

  EXPECT_EQ(std::filesystem::last_write_time(path), before_time);
  EXPECT_EQ(dret::DiffNoTreeCacheKey(config_.keys, dret::kDiffBlockSize), key);
  // ...and no second copy under a block-size key.
  EXPECT_FALSE(std::filesystem::exists(
      sdsl::cache_file_name<TDiff>(std::format("bs{}_{}", dret::kDiffBlockSize, key), config_)));
}

// A PDL-RP core cached with the pre-2026-10 codec (RPCodecWithLengths, whose
// lists carry a per-rule length array they never read) is converted on
// construction instead of rebuilt: the new core answers the same and is no
// larger, and the legacy file is left alone.
TEST_F(RMQSLPCacheReuseTest, pdl_rp_converts_legacy_core) {
  using S = ExternalGenericStorage;
  using Legacy = dret::pdl::DocListIdxPDL<S, dret::Alphabet<>, sri::RIndexCount<S, dret::Alphabet<>>,
                                          dret::pdl::PDLGetDocsDA<S>, dret::pdl::RPCodecWithLengths>;
  using Current = dret::pdl::DocListIdxPDLRP<S>;
  const std::vector<std::pair<std::string, std::vector<std::size_t>>> expected = {
      {"TAT", {0}}, {"LAT", {1}}, {"LAL", {2}}, {"TA", {0, 1}},
      {"LA", {1, 2}}, {"A", {0, 1, 2}}, {"TAL", {}},
  };

  Legacy legacy(std::ref(storage_), 2, 2.0f);
  construct(legacy, config_);

  Current current(std::ref(storage_), 2, 2.0f);
  construct(current, config_);

  EXPECT_LE(sdsl::size_in_bytes(current), sdsl::size_in_bytes(legacy));
  for (const auto& [pattern, docs] : expected) {
    DocListResultVector result;
    current.Search(pattern, std::ref(result));
    result();
    EXPECT_THAT(result, testing::ElementsAreArray(docs)) << pattern;
  }
}

// GCDA-nolists' differential index sweeps the same sample spacing as the RMQ
// backend: each spacing answers correctly, keeps its own file (512 the one
// GCDA-nolists always used), and all of them share one base grammar.
TEST_F(RMQSLPCacheReuseTest, gcda_nolists_diff_sweeps_spacing) {
  using S = ExternalGenericStorage;
  using TDiff = dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>;
  using Index = dret::DocListIdxSLP<S, dret::Alphabet<>, RMQCountIdx<S>, TDiff>;
  const std::vector<std::pair<std::string, std::vector<std::size_t>>> expected = {
      {"TAT", {0}}, {"LAT", {1}}, {"LAL", {2}}, {"TA", {0, 1}},
      {"LA", {1, 2}}, {"A", {0, 1, 2}}, {"TAL", {}},
  };
  std::vector<std::size_t> sizes;
  for (std::uint32_t spacing : {2u, 4u, 512u}) {
    Index index(std::ref(storage_), spacing);
    construct(index, config_);
    sizes.push_back(sdsl::size_in_bytes(index));
    for (const auto& [pattern, docs] : expected) {
      DocListResultVector result;
      index.Search(pattern, std::ref(result));
      result();
      EXPECT_THAT(result, testing::ElementsAreArray(docs)) << "spacing " << spacing << " " << pattern;
    }
    const auto key = dret::DiffNoTreeCacheKey(config_.keys, spacing);
    EXPECT_EQ(key == config_.keys[dret::conf::kDSLPNS].get<std::string>(), spacing == dret::kDiffBlockSize);
    EXPECT_TRUE(std::filesystem::exists(sdsl::cache_file_name<TDiff>(key, config_))) << key;
  }
  const auto dir = std::filesystem::path(
      sdsl::cache_file_name<TDiff>(dret::DiffNoTreeCacheKey(config_.keys, 2u), config_)).parent_path();
  std::size_t n_grammars = 0;
  for (const auto& e : std::filesystem::directory_iterator(dir))
    n_grammars += e.path().filename().string().find(dret::DiffNoTreeGrammarKey(config_.keys)) != std::string::npos;
  EXPECT_EQ(n_grammars, 1u);
  // The spacing is real: a denser sampling stores more samples.
  EXPECT_GT(sizes[0], sizes[2]);
}

// PrefixedKey reproduces the parameterised file names dret used to format by
// hand, so moving a call site to it renames nothing.
TEST(PrefixedKeyTest, matches_the_formats_in_use) {
  const auto keys = dret::Keys<8>().keys;
  EXPECT_EQ(dret::PrefixedKey(keys, dret::conf::kBsSf, "gcda_slp", 512u, 4.0f), std::format("{}-{}_gcda_slp", 512u, 4.0f));
  EXPECT_EQ(dret::PrefixedKey(keys, dret::conf::kBsSf, "gcda_docs", 1024u, 32.0f), "1024-32_gcda_docs");
  EXPECT_EQ(dret::PrefixedKey(keys, dret::conf::kSpacing, "dslp_ns", 128u), "bs128_dslp_ns");
  EXPECT_EQ(dret::DiffNoTreeCacheKey(keys, 128u), "bs128_dslp_ns");
  EXPECT_EQ(dret::DiffNoTreeCacheKey(keys, 512u), "dslp_ns");
  EXPECT_EQ(dret::DiffNoTreeGrammarKey(keys), "dslp_ns_grammar");
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
