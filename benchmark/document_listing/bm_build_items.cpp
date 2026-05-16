//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/8/21.
//

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/vlc_vector.hpp>

#include "../bm_size_counters.h"

#include "dret/slp/basic_slp_span_length.h"
#include "dret/construct_base.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/doc_list/doc_list_slp.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/doc_list/doc_list_dgcda.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/index_base.h"
#include "dret/doc_list/doc_list_pdl_bc.h"
#include "dret/doc_list/doc_list_pdl_plain.h"
#include "dret/doc_list/doc_list_pdl_rp.h"
#include "dret/pdl/get_docs.h"
#include "dret/pdl/storage_policy.h"
#include "dret/size_report.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_int32(data_width, 8, "Data width in bits: 8, 16, 32 or 64");
DEFINE_string(sa_algo, "SDSL_SE_SAIS", "Suffix Array Algorithm: SDSL_SE_SAIS, SDSL_LIBDIVSUFSORT, BIG_BWT");
DEFINE_int32(doc_delim, 0, "Document delimiter.");

DEFINE_int32(min_block_size, 512, "Minimum block size (power of 2).");
DEFINE_int32(max_block_size, 512, "Maximum block size (power of 2).");
DEFINE_int32(min_storing_factor, 4, "Minimum storing factor (power of 2).");
DEFINE_int32(max_storing_factor, 4, "Maximum storing factor (power of 2).");
DEFINE_string(rmq_get_doc_variants, "da,slp,slp_ns", "RMQ GetDoc variants to build: comma-separated da,slp,dslp.");

DEFINE_string(gcda_slp_variants,
              "default,compact_bp,compact_louds,cslp",
              "GCDA TSLP variants: comma-separated default,compact_bp,compact_louds,cslp.");

DEFINE_string(bare_slp_variants,
              "default,raw,dv,vv",
              "Bare-SLP container variants for SLP-NS family: comma-separated default,raw,dv,vv.");

DEFINE_string(pdl_variants,
              "",
              "PDL stored-set codec variants: comma-separated plain,rp,bc. Empty disables PDL.");
DEFINE_string(pdl_get_doc_variants,
              "da",
              "PDL raw-range get-doc variants: comma-separated da,slp,dslp.");
DEFINE_string(pdl_storage_policy,
              "occurrence_weighted",
              "PDL storage policy: comma-separated "
              "occurrence_weighted,all_internal,leaves_only.");

//~~~~~~~

enum class RMQGetDocVariant {
  DA,
  SLP,
  SLP_NS,
  DSLP,
};

std::vector<RMQGetDocVariant> ParseRMQGetDocVariants(const std::string& value) {
  std::vector<RMQGetDocVariant> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "da") {
      variants.push_back(RMQGetDocVariant::DA);
    } else if (item == "slp") {
      variants.push_back(RMQGetDocVariant::SLP);
    } else if (item == "slp_ns") {
      variants.push_back(RMQGetDocVariant::SLP_NS);
    } else if (item == "dslp") {
      variants.push_back(RMQGetDocVariant::DSLP);
    } else if (!item.empty()) {
      throw std::invalid_argument("Unknown --rmq_get_doc_variants item: " + item);
    }
  }
  if (variants.empty())
    variants.push_back(RMQGetDocVariant::DA);
  return variants;
}

enum class GCDASLPVariant {
  Default,
  CompactBP,
  CompactLOUDS,
  CSLP,
};

enum class BareSLPVariant {
  Default,
  Raw,
  DV,
  VV,
};

std::vector<BareSLPVariant> ParseBareSLPVariants(const std::string& value) {
  std::vector<BareSLPVariant> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "default") {
      variants.push_back(BareSLPVariant::Default);
    } else if (item == "raw") {
      variants.push_back(BareSLPVariant::Raw);
    } else if (item == "dv") {
      variants.push_back(BareSLPVariant::DV);
    } else if (item == "vv") {
      variants.push_back(BareSLPVariant::VV);
    } else if (!item.empty()) {
      throw std::invalid_argument("Unknown --bare_slp_variants item: " + item);
    }
  }
  if (variants.empty())
    variants.push_back(BareSLPVariant::Default);
  return variants;
}

std::vector<GCDASLPVariant> ParseGCDASLPVariants(const std::string& value) {
  std::vector<GCDASLPVariant> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "default") {
      variants.push_back(GCDASLPVariant::Default);
    } else if (item == "compact_bp") {
      variants.push_back(GCDASLPVariant::CompactBP);
    } else if (item == "compact_louds") {
      variants.push_back(GCDASLPVariant::CompactLOUDS);
    } else if (item == "cslp") {
      variants.push_back(GCDASLPVariant::CSLP);
    } else if (!item.empty()) {
      throw std::invalid_argument("Unknown --gcda_slp_variants item: " + item);
    }
  }
  if (variants.empty())
    variants.push_back(GCDASLPVariant::Default);
  return variants;
}

bool HasVariant(const std::vector<GCDASLPVariant>& variants, GCDASLPVariant variant) {
  return std::find(variants.begin(), variants.end(), variant) != variants.end();
}

bool HasVariant(const std::vector<RMQGetDocVariant>& variants, RMQGetDocVariant variant) {
  return std::find(variants.begin(), variants.end(), variant) != variants.end();
}

bool HasVariant(const std::vector<BareSLPVariant>& variants, BareSLPVariant variant) {
  return std::find(variants.begin(), variants.end(), variant) != variants.end();
}

enum class PDLCodecVariant { Plain, RP, BC };
enum class PDLGetDocVariant { DA, SLP, DSLP };

std::vector<PDLCodecVariant> ParsePDLVariants(const std::string& value) {
  std::vector<PDLCodecVariant> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "plain")
      variants.push_back(PDLCodecVariant::Plain);
    else if (item == "rp")
      variants.push_back(PDLCodecVariant::RP);
    else if (item == "bc")
      variants.push_back(PDLCodecVariant::BC);
    else if (!item.empty())
      throw std::invalid_argument("Unknown --pdl_variants item: " + item);
  }
  return variants;
}

std::vector<PDLGetDocVariant> ParsePDLGetDocVariants(const std::string& value) {
  std::vector<PDLGetDocVariant> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "da")
      variants.push_back(PDLGetDocVariant::DA);
    else if (item == "slp")
      variants.push_back(PDLGetDocVariant::SLP);
    else if (item == "dslp")
      variants.push_back(PDLGetDocVariant::DSLP);
    else if (!item.empty())
      throw std::invalid_argument("Unknown --pdl_get_doc_variants item: " + item);
  }
  if (variants.empty())
    variants.push_back(PDLGetDocVariant::DA);
  return variants;
}

std::vector<dret::pdl::StoragePolicy> ParsePDLStoragePolicies(const std::string& value) {
  std::vector<dret::pdl::StoragePolicy> policies;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "occurrence_weighted")
      policies.push_back(dret::pdl::StoragePolicy::OccurrenceWeighted);
    else if (item == "all_internal")
      policies.push_back(dret::pdl::StoragePolicy::StoreAllInternal);
    else if (item == "leaves_only")
      policies.push_back(dret::pdl::StoragePolicy::LeavesOnly);
    else if (!item.empty())
      throw std::invalid_argument("Unknown --pdl_storage_policy item: " + item);
  }
  if (policies.empty())
    policies.push_back(dret::pdl::StoragePolicy::OccurrenceWeighted);
  return policies;
}

const char* PDLCodecVariantName(PDLCodecVariant v) {
  switch (v) {
    case PDLCodecVariant::Plain: return "Plain";
    case PDLCodecVariant::RP:    return "RP";
    case PDLCodecVariant::BC:    return "BC";
  }
  return "UNKNOWN";
}

const char* PDLGetDocVariantName(PDLGetDocVariant v) {
  switch (v) {
    case PDLGetDocVariant::DA:   return "DA";
    case PDLGetDocVariant::SLP:  return "SLP";
    case PDLGetDocVariant::DSLP: return "DSLP";
  }
  return "UNKNOWN";
}

const char* PDLStoragePolicyName(dret::pdl::StoragePolicy p) {
  switch (p) {
    case dret::pdl::StoragePolicy::OccurrenceWeighted: return "OccurrenceWeighted";
    case dret::pdl::StoragePolicy::StoreAllInternal:   return "StoreAllInternal";
    case dret::pdl::StoragePolicy::LeavesOnly:         return "LeavesOnly";
  }
  return "UNKNOWN";
}

//~~~~~~~


void SetupCommonCounters(benchmark::State& t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["bs"] = 0;  // Block size
  t_state.counters["sf"] = 0;  // Storing factor
}

//~~~~~~~


template <typename TIndex>
void BM_ConstructBruteIdx(benchmark::State& t_state, dret::Config t_config, const std::string& t_data_path) {
  TIndex index;

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }

  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  {
    std::ofstream ofs("construction-" + idx_name + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  }
  {
    std::ofstream ofs("construction-" + idx_name + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  }

  SetupCommonCounters(t_state);
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }
}

//~~~~~~~


template <typename TIndex>
void BM_ConstructDocListIdxGCDA(benchmark::State& t_state, dret::Config t_config) {
  uint32_t block_size = static_cast<uint32_t>(t_state.range(0));
  float storing_factor = static_cast<float>(t_state.range(1));

  dret::GenericStorage storage;
  TIndex index(storage, block_size, storing_factor);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }

  auto suffix = "-bs" + std::to_string(block_size) + "-sf" + std::to_string(static_cast<int>(storing_factor));
  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(0, bm_name.find('/'));
  {
    std::ofstream ofs("construction-" + idx_name + suffix + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  }
  {
    std::ofstream ofs("construction-" + idx_name + suffix + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  }

  SetupCommonCounters(t_state);
  t_state.counters["bs"] = block_size;
  t_state.counters["sf"] = storing_factor;
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }

  // Per-field size breakdown: populates benchmark counters and writes a JSON sidecar.
  // Built index must be loaded (not just cache-populated) so GetSizeReport can dispatch.
  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  dret::writeSizesJson("sizes-" + idx_name + suffix + ".json", sizes);
}

//~~~~~~~

// Phase C: parameter-less GCDA-shaped construct benchmark for the
// non-sampled `DocListIdxSLP`. No block_size / storing_factor knobs.
template <typename TIndex>
void BM_ConstructDocListIdxSLP(benchmark::State& t_state, dret::Config t_config) {
  dret::GenericStorage storage;
  TIndex index(storage);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }

  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  {
    std::ofstream ofs("construction-" + idx_name + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  }
  {
    std::ofstream ofs("construction-" + idx_name + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  }

  SetupCommonCounters(t_state);
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }

  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  dret::writeSizesJson("sizes-" + idx_name + ".json", sizes);
}

//~~~~~~~

// Parameter-less RMQ-Compressed harness for SLP-NS variants. The bare grammar::SLP<>
// is parameter-free (no bs/sf knobs), so the RMQ core is constructed single-arg.
template <typename TIndex, typename TCore>
void BM_ConstructDocListIdxRMQCompressedNS(benchmark::State& t_state, dret::Config t_config) {
  dret::GenericStorage storage;
  TCore core(storage);
  TIndex index(storage, core);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }

  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  {
    std::ofstream ofs("construction-" + idx_name + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  }
  {
    std::ofstream ofs("construction-" + idx_name + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  }

  SetupCommonCounters(t_state);
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }

  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  dret::writeSizesJson("sizes-" + idx_name + ".json", sizes);
}

//~~~~~~~

template <typename TIndex, typename TCore>
void BM_ConstructDocListIdxRMQCompressed(benchmark::State& t_state, dret::Config t_config) {
  uint32_t block_size = static_cast<uint32_t>(t_state.range(0));
  float storing_factor = static_cast<float>(t_state.range(1));

  dret::GenericStorage storage;
  TCore core(storage, block_size, storing_factor);
  TIndex index(storage, core);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }

  SetupCommonCounters(t_state);
  t_state.counters["bs"] = block_size;
  t_state.counters["sf"] = storing_factor;
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }

  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
}

//~~~~~~~


// PDL takes a storage_policy argument in addition to block_size /
// storing_factor, and is parameterized on a (codec, get-doc) pair at
// the type level. This template instantiates one such pair.
template <typename TIndex>
void BM_ConstructDocListIdxPDL(benchmark::State& t_state,
                                dret::Config t_config,
                                dret::pdl::StoragePolicy t_policy) {
  uint32_t block_size = static_cast<uint32_t>(t_state.range(0));
  float storing_factor = static_cast<float>(t_state.range(1));

  dret::GenericStorage storage;
  TIndex index(storage, block_size, storing_factor, t_policy);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }

  auto suffix = "-bs" + std::to_string(block_size) + "-sf" + std::to_string(static_cast<int>(storing_factor));
  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(0, bm_name.find('/'));
  {
    std::ofstream ofs("construction-" + idx_name + suffix + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  }
  {
    std::ofstream ofs("construction-" + idx_name + suffix + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  }

  SetupCommonCounters(t_state);
  t_state.counters["bs"] = block_size;
  t_state.counters["sf"] = storing_factor;
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }

  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  dret::writeSizesJson("sizes-" + idx_name + suffix + ".json", sizes);
}

//~~~~~~~


static std::vector<int64_t> powersOfTwo(int32_t t_min, int32_t t_max) {
  std::vector<int64_t> result;
  for (int64_t v = t_min; v <= t_max; v *= 2)
    result.push_back(v);
  return result;
}

//~~~~~~~


int main(int argc, char** argv) {
  gflags::SetUsageMessage("This program calculates the ri items for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  std::string data_path = FLAGS_data;
  std::vector<RMQGetDocVariant> rmq_get_doc_variants;
  std::vector<GCDASLPVariant> gcda_slp_variants;
  std::vector<BareSLPVariant> bare_slp_variants;
  std::vector<PDLCodecVariant> pdl_variants;
  std::vector<PDLGetDocVariant> pdl_get_doc_variants;
  std::vector<dret::pdl::StoragePolicy> pdl_storage_policies;
  try {
    rmq_get_doc_variants = ParseRMQGetDocVariants(FLAGS_rmq_get_doc_variants);
    gcda_slp_variants = ParseGCDASLPVariants(FLAGS_gcda_slp_variants);
    bare_slp_variants = ParseBareSLPVariants(FLAGS_bare_slp_variants);
    pdl_variants = ParsePDLVariants(FLAGS_pdl_variants);
    pdl_get_doc_variants = ParsePDLGetDocVariants(FLAGS_pdl_get_doc_variants);
    pdl_storage_policies = ParsePDLStoragePolicies(FLAGS_pdl_storage_policy);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  dret::Config config(data_path,
                      std::filesystem::current_path(),
                      sri::toSAAlgo(FLAGS_sa_algo),
                      false,
                      FLAGS_data_width,
                      FLAGS_doc_delim);

  benchmark::RegisterBenchmark("DocListIdxBrute", BM_ConstructBruteIdx<dret::DocListIdxBrute<>>, config, data_path);

  using SADAIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                            dret::Alphabet<>,
                                            sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                            dret::rmq::SadaCore<dret::GenericStorage>>;
  using ILCPIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                            dret::Alphabet<>,
                                            sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                            dret::rmq::IlcpCore<dret::GenericStorage>>;
  using CILCPIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                             dret::Alphabet<>,
                                             sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                             dret::rmq::CilcpCore<dret::GenericStorage>>;
  if (HasVariant(rmq_get_doc_variants, RMQGetDocVariant::DA)) {
    benchmark::RegisterBenchmark("DocListSADA-DA", BM_ConstructBruteIdx<SADAIdx>, config, data_path);
    benchmark::RegisterBenchmark("DocListILCP-DA", BM_ConstructBruteIdx<ILCPIdx>, config, data_path);
    benchmark::RegisterBenchmark("DocListCILCP-DA", BM_ConstructBruteIdx<CILCPIdx>, config, data_path);
  }

  // Phase C: non-sampled SLP index, plus the three RMQ listing cores backed by
  // the bare grammar::SLP<> cache (kSLPNS). All four are parameter-free (no bs/sf
  // axis). Each gets fanned out across the bare-SLP container axis (Default =
  // sdsl::int_vector<>, DV = sdsl::dac_vector<>, VV = sdsl::vlc_vector<>); type-
  // hashes on grammar::SLP<TVars,TLens> distinguish the on-disk cache files.
  auto register_slp_ns_for_tslp = [&]<typename TSLP>(const char* suffix) {
    using IdxSLP = dret::DocListIdxSLP<dret::GenericStorage,
                                        dret::Alphabet<>,
                                        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                        TSLP>;
    using GetDoc = dret::rmq::GetDocSLP_NS<dret::GenericStorage,
                                            dret::Alphabet<>::int_width,
                                            TSLP>;
    using SadaCore = dret::rmq::SadaCore<dret::GenericStorage,
                                          dret::Alphabet<>::int_width,
                                          sdsl::rmq_succinct_sct<true>,
                                          sdsl::sd_vector<>,
                                          GetDoc>;
    using IlcpCore = dret::rmq::IlcpCore<dret::GenericStorage,
                                          dret::Alphabet<>::int_width,
                                          sdsl::sd_vector<>,
                                          sdsl::rmq_succinct_sct<true>,
                                          sdsl::sd_vector<>,
                                          GetDoc>;
    using CilcpCore = dret::rmq::CilcpCore<dret::GenericStorage,
                                            dret::Alphabet<>::int_width,
                                            sdsl::sd_vector<>,
                                            sdsl::rmq_succinct_sct<true>,
                                            sdsl::sd_vector<>,
                                            GetDoc>;
    using SadaIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                              dret::Alphabet<>,
                                              sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                              SadaCore>;
    using IlcpIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                              dret::Alphabet<>,
                                              sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                              IlcpCore>;
    using CilcpIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                               dret::Alphabet<>,
                                               sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                               CilcpCore>;

    benchmark::RegisterBenchmark(
        std::string("DocListSLP-NS") + suffix, BM_ConstructDocListIdxSLP<IdxSLP>, config);
    if (HasVariant(rmq_get_doc_variants, RMQGetDocVariant::SLP_NS)) {
      benchmark::RegisterBenchmark(std::string("DocListSADA-SLP-NS") + suffix,
                                   BM_ConstructDocListIdxRMQCompressedNS<SadaIdx, SadaCore>, config);
      benchmark::RegisterBenchmark(std::string("DocListILCP-SLP-NS") + suffix,
                                   BM_ConstructDocListIdxRMQCompressedNS<IlcpIdx, IlcpCore>, config);
      benchmark::RegisterBenchmark(std::string("DocListCILCP-SLP-NS") + suffix,
                                   BM_ConstructDocListIdxRMQCompressedNS<CilcpIdx, CilcpCore>, config);
    }
  };

  if (HasVariant(bare_slp_variants, BareSLPVariant::Default))
    register_slp_ns_for_tslp.template operator()<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>("");
  if (HasVariant(bare_slp_variants, BareSLPVariant::Raw))
    register_slp_ns_for_tslp.template operator()<grammar::SLP<>>("-Raw");
  if (HasVariant(bare_slp_variants, BareSLPVariant::DV))
    register_slp_ns_for_tslp.template operator()<grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>>("-DV");
  if (HasVariant(bare_slp_variants, BareSLPVariant::VV))
    register_slp_ns_for_tslp.template operator()<grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>>("-VV");

  auto block_sizes = powersOfTwo(FLAGS_min_block_size, FLAGS_max_block_size);
  auto storing_factors = powersOfTwo(FLAGS_min_storing_factor, FLAGS_max_storing_factor);
  if (!block_sizes.empty() && !storing_factors.empty()) {
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::Default)) {
      benchmark::RegisterBenchmark(
          "DocListGCDA", BM_ConstructDocListIdxGCDA<dret::gcda::DocListIdxGCDA<>>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    }

    using GCDA_CompactBP = dret::gcda::DocListIdxGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        grammar::CompactBPSLP<>>;
    using GCDA_CompactLOUDS = dret::gcda::DocListIdxGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        grammar::CompactLOUDSSLP<>>;
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactBP)) {
      benchmark::RegisterBenchmark(
          "DocListGCDA-CompactBP", BM_ConstructDocListIdxGCDA<GCDA_CompactBP>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    }
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactLOUDS)) {
      benchmark::RegisterBenchmark(
          "DocListGCDA-CompactLOUDS", BM_ConstructDocListIdxGCDA<GCDA_CompactLOUDS>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    }

    using GCDA_CSLP = dret::gcda::DocListIdxGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        grammar::CombinedSLPWithUnitCover<>>;
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::CSLP)) {
      benchmark::RegisterBenchmark(
          "DocListGCDA-CSLP", BM_ConstructDocListIdxGCDA<GCDA_CSLP>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    }

    benchmark::RegisterBenchmark(
        "DocListDGCDA", BM_ConstructDocListIdxGCDA<dret::dgcda::DocListIdxDGCDA<>>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    using GetDocSLP = dret::rmq::GetDocSLP<dret::GenericStorage>;
    using SADACoreSLP = dret::rmq::SadaCore<dret::GenericStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             GetDocSLP>;
    using ILCPCoreSLP = dret::rmq::IlcpCore<dret::GenericStorage,
                                             dret::Alphabet<>::int_width,
                                             sdsl::sd_vector<>,
                                             sdsl::rmq_succinct_sct<true>,
                                             sdsl::sd_vector<>,
                                             GetDocSLP>;
    using CILCPCoreSLP = dret::rmq::CilcpCore<dret::GenericStorage,
                                               dret::Alphabet<>::int_width,
                                               sdsl::sd_vector<>,
                                               sdsl::rmq_succinct_sct<true>,
                                               sdsl::sd_vector<>,
                                               GetDocSLP>;
    using SADAIdxSLP = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                dret::Alphabet<>,
                                                sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                SADACoreSLP>;
    using ILCPIdxSLP = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                dret::Alphabet<>,
                                                sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                ILCPCoreSLP>;
    using CILCPIdxSLP = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                 dret::Alphabet<>,
                                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                 CILCPCoreSLP>;
    if (HasVariant(rmq_get_doc_variants, RMQGetDocVariant::SLP)) {
      // Default SLP TSLP (matches GCDA default).
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::Default)) {
        benchmark::RegisterBenchmark("DocListSADA-SLP", BM_ConstructDocListIdxRMQCompressed<SADAIdxSLP, SADACoreSLP>, config)
            ->ArgsProduct({block_sizes, storing_factors});
        benchmark::RegisterBenchmark("DocListILCP-SLP", BM_ConstructDocListIdxRMQCompressed<ILCPIdxSLP, ILCPCoreSLP>, config)
            ->ArgsProduct({block_sizes, storing_factors});
        benchmark::RegisterBenchmark("DocListCILCP-SLP", BM_ConstructDocListIdxRMQCompressed<CILCPIdxSLP, CILCPCoreSLP>, config)
            ->ArgsProduct({block_sizes, storing_factors});
      }

      // Compact-grammar TSLPs — share the SLP cache file with the matching
      // DocListGCDA-Compact* build via TSLP type matching.
      auto register_rmq_slp_for_tslp = [&]<typename TSLP>(const char* suffix) {
        using GetDocSLPVar = dret::rmq::GetDocSLP<dret::GenericStorage,
                                                  dret::Alphabet<>::int_width,
                                                  TSLP>;
        using SadaCore = dret::rmq::SadaCore<dret::GenericStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              GetDocSLPVar>;
        using IlcpCore = dret::rmq::IlcpCore<dret::GenericStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::sd_vector<>,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              GetDocSLPVar>;
        using CilcpCore = dret::rmq::CilcpCore<dret::GenericStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                GetDocSLPVar>;
        using SadaIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                  dret::Alphabet<>,
                                                  sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                  SadaCore>;
        using IlcpIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                  dret::Alphabet<>,
                                                  sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                  IlcpCore>;
        using CilcpIdx = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                   dret::Alphabet<>,
                                                   sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                   CilcpCore>;
        benchmark::RegisterBenchmark(std::string("DocListSADA-SLP-") + suffix,
                                     BM_ConstructDocListIdxRMQCompressed<SadaIdx, SadaCore>, config)
            ->ArgsProduct({block_sizes, storing_factors});
        benchmark::RegisterBenchmark(std::string("DocListILCP-SLP-") + suffix,
                                     BM_ConstructDocListIdxRMQCompressed<IlcpIdx, IlcpCore>, config)
            ->ArgsProduct({block_sizes, storing_factors});
        benchmark::RegisterBenchmark(std::string("DocListCILCP-SLP-") + suffix,
                                     BM_ConstructDocListIdxRMQCompressed<CilcpIdx, CilcpCore>, config)
            ->ArgsProduct({block_sizes, storing_factors});
      };

      if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactBP))
        register_rmq_slp_for_tslp.template operator()<grammar::CompactBPSLP<>>("CompactBP");
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactLOUDS))
        register_rmq_slp_for_tslp.template operator()<grammar::CompactLOUDSSLP<>>("CompactLOUDS");
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::CSLP))
        register_rmq_slp_for_tslp.template operator()<grammar::CombinedSLPWithUnitCover<>>("CSLP");
    }

    using GetDocDSLP = dret::rmq::GetDocDSLP<dret::GenericStorage>;
    using SADACoreDSLP = dret::rmq::SadaCore<dret::GenericStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              GetDocDSLP>;
    using ILCPCoreDSLP = dret::rmq::IlcpCore<dret::GenericStorage,
                                              dret::Alphabet<>::int_width,
                                              sdsl::sd_vector<>,
                                              sdsl::rmq_succinct_sct<true>,
                                              sdsl::sd_vector<>,
                                              GetDocDSLP>;
    using CILCPCoreDSLP = dret::rmq::CilcpCore<dret::GenericStorage,
                                                dret::Alphabet<>::int_width,
                                                sdsl::sd_vector<>,
                                                sdsl::rmq_succinct_sct<true>,
                                                sdsl::sd_vector<>,
                                                GetDocDSLP>;
    using SADAIdxDSLP = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                 dret::Alphabet<>,
                                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                 SADACoreDSLP>;
    using ILCPIdxDSLP = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                 dret::Alphabet<>,
                                                 sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                 ILCPCoreDSLP>;
    using CILCPIdxDSLP = dret::rmq::DocListIdxRMQ<dret::GenericStorage,
                                                  dret::Alphabet<>,
                                                  sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
                                                  CILCPCoreDSLP>;
    if (HasVariant(rmq_get_doc_variants, RMQGetDocVariant::DSLP)) {
      benchmark::RegisterBenchmark("DocListSADA-DSLP", BM_ConstructDocListIdxRMQCompressed<SADAIdxDSLP, SADACoreDSLP>, config)
          ->ArgsProduct({block_sizes, storing_factors});
      benchmark::RegisterBenchmark("DocListILCP-DSLP", BM_ConstructDocListIdxRMQCompressed<ILCPIdxDSLP, ILCPCoreDSLP>, config)
          ->ArgsProduct({block_sizes, storing_factors});
      benchmark::RegisterBenchmark("DocListCILCP-DSLP", BM_ConstructDocListIdxRMQCompressed<CILCPIdxDSLP, CILCPCoreDSLP>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    }

    using DGCDA_OTF = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-OTF", BM_ConstructDocListIdxGCDA<DGCDA_OTF>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    using DGCDA_CRL = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-CRL", BM_ConstructDocListIdxGCDA<DGCDA_CRL>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    // Compressed int-vector variants: vary TRoots/TSpanSums/TSamples in lockstep;
    // TSampleRootsPos stays at the class default sdsl::enc_vector<>.
    using DGCDA_EV = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        dret::DifferentialLightSLP<grammar::SLP<>,
                                    grammar::SampledSLP<>,
                                    sdsl::enc_vector<>,
                                    sdsl::enc_vector<>,
                                    sdsl::enc_vector<>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-EV", BM_ConstructDocListIdxGCDA<DGCDA_EV>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    using DGCDA_DV = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        dret::DifferentialLightSLP<grammar::SLP<>,
                                    grammar::SampledSLP<>,
                                    sdsl::dac_vector<>,
                                    sdsl::dac_vector<>,
                                    sdsl::dac_vector<>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-DV", BM_ConstructDocListIdxGCDA<DGCDA_DV>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    using DGCDA_VV = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>,
        dret::DifferentialLightSLP<grammar::SLP<>,
                                    grammar::SampledSLP<>,
                                    sdsl::vlc_vector<>,
                                    sdsl::vlc_vector<>,
                                    sdsl::vlc_vector<>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-VV", BM_ConstructDocListIdxGCDA<DGCDA_VV>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    // PDL build benchmarks. The (codec, get-doc) pair is a type-level
    // axis (9 template instantiations); storage policy is a runtime
    // arg threaded into the index constructor. Skipped entirely when
    // --pdl_variants is empty.
    auto register_pdl = [&]<typename TIndex>(const std::string& base) {
      for (const auto sp : pdl_storage_policies) {
        std::string name = "DocListPDL-" + base + "-" + PDLStoragePolicyName(sp);
        benchmark::RegisterBenchmark(
            name, BM_ConstructDocListIdxPDL<TIndex>, config, sp)
            ->ArgsProduct({block_sizes, storing_factors});
      }
    };

    for (const auto pdl_v : pdl_variants) {
      for (const auto pdl_gd : pdl_get_doc_variants) {
        std::string pair = std::string(PDLCodecVariantName(pdl_v)) + "-" +
                           std::string(PDLGetDocVariantName(pdl_gd));
        using TCountIdx = sri::RIndexCount<dret::GenericStorage, dret::Alphabet<>>;
        switch (pdl_v) {
          case PDLCodecVariant::Plain:
            switch (pdl_gd) {
              case PDLGetDocVariant::DA:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLPlain<dret::GenericStorage>>(pair);
                break;
              case PDLGetDocVariant::SLP:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLPlain<dret::GenericStorage, dret::Alphabet<>, TCountIdx,
                                                  dret::pdl::PDLGetDocsSLP<dret::GenericStorage>>>(pair);
                break;
              case PDLGetDocVariant::DSLP:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLPlain<dret::GenericStorage, dret::Alphabet<>, TCountIdx,
                                                  dret::pdl::PDLGetDocsDSLP<dret::GenericStorage>>>(pair);
                break;
            }
            break;
          case PDLCodecVariant::RP:
            switch (pdl_gd) {
              case PDLGetDocVariant::DA:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLRP<dret::GenericStorage>>(pair);
                break;
              case PDLGetDocVariant::SLP:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLRP<dret::GenericStorage, dret::Alphabet<>, TCountIdx,
                                               dret::pdl::PDLGetDocsSLP<dret::GenericStorage>>>(pair);
                break;
              case PDLGetDocVariant::DSLP:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLRP<dret::GenericStorage, dret::Alphabet<>, TCountIdx,
                                               dret::pdl::PDLGetDocsDSLP<dret::GenericStorage>>>(pair);
                break;
            }
            break;
          case PDLCodecVariant::BC:
            switch (pdl_gd) {
              case PDLGetDocVariant::DA:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLBC<dret::GenericStorage>>(pair);
                break;
              case PDLGetDocVariant::SLP:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLBC<dret::GenericStorage, dret::Alphabet<>, TCountIdx,
                                               dret::pdl::PDLGetDocsSLP<dret::GenericStorage>>>(pair);
                break;
              case PDLGetDocVariant::DSLP:
                register_pdl.template operator()<
                    dret::pdl::DocListIdxPDLBC<dret::GenericStorage, dret::Alphabet<>, TCountIdx,
                                               dret::pdl::PDLGetDocsDSLP<dret::GenericStorage>>>(pair);
                break;
            }
            break;
        }
      }
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
