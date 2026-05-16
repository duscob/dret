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
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/index_base.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"
#include "dret/pdl/storage_policy.h"
#include "dret/size_report.h"

#include "enum_traits.h"
#include "factories/brute.h"
#include "factories/dgcda.h"
#include "factories/gcda.h"
#include "factories/pdl.h"
#include "factories/rmq.h"
#include "factories/slp_ns.h"


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
              "light,compact_bp,compact_louds,combined",
              "GCDA TSLP variants: comma-separated light,compact_bp,compact_louds,combined.");

DEFINE_string(bare_slp_variants,
              "iv,raw,dv,vv",
              "Bare-SLP container variants for SLP-NS family: comma-separated iv,raw,dv,vv.");

DEFINE_string(dgcda_slp_variants,
              "default,otf,crl,ev,dv,vv",
              "DGCDA TSLP variants: comma-separated default,otf,crl,ev,dv,vv.");

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

// All axis enums and parsing/name helpers come from enum_traits.h. The local
// HasVariant<T> alias preserves the call sites further below.
using bench::axes::BareSLPVariant;
using bench::axes::DGCDASLPVariant;
using bench::axes::GCDASLPVariant;
using bench::axes::GetDocEnum;
using bench::axes::PDLStoragePolicy;
using bench::axes::PDLVariant;

template <typename E>
bool HasVariant(const std::vector<E>& vec, E v) {
  return bench::axes::Contains(vec, v);
}

// PDL accepts only DA/SLP/DSLP (no SLP_NS). The generic parser allows all
// four GetDocEnum values, so we filter SLP_NS afterwards here.
static std::vector<GetDocEnum> ParseGetDocEnumsStrict(const std::string& value) {
  auto variants = bench::axes::ParseCSV<GetDocEnum>(value);
  for (auto v : variants) {
    if (v == GetDocEnum::SLP_NS) {
      throw std::invalid_argument("--pdl_get_doc_variants: 'slp_ns' is not a valid PDL backing");
    }
  }
  return variants;
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
  std::vector<GetDocEnum> rmq_get_doc_variants;
  std::vector<GCDASLPVariant> gcda_slp_variants;
  std::vector<BareSLPVariant> bare_slp_variants;
  std::vector<DGCDASLPVariant> dgcda_slp_variants;
  std::vector<PDLVariant> pdl_variants;
  std::vector<GetDocEnum> pdl_get_doc_variants;
  std::vector<PDLStoragePolicy> pdl_storage_policies;
  try {
    rmq_get_doc_variants = bench::axes::ParseCSV<GetDocEnum>(FLAGS_rmq_get_doc_variants);
    gcda_slp_variants = bench::axes::ParseCSV<GCDASLPVariant>(FLAGS_gcda_slp_variants);
    bare_slp_variants = bench::axes::ParseCSV<BareSLPVariant>(FLAGS_bare_slp_variants);
    dgcda_slp_variants = bench::axes::ParseCSV<DGCDASLPVariant>(FLAGS_dgcda_slp_variants);
    pdl_variants = bench::axes::ParseCSV<PDLVariant>(FLAGS_pdl_variants);
    pdl_get_doc_variants = ParseGetDocEnumsStrict(FLAGS_pdl_get_doc_variants);
    pdl_storage_policies = bench::axes::ParseCSV<PDLStoragePolicy>(FLAGS_pdl_storage_policy);
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

  // Brute baseline — DocListIdxBrute<> with default storage. (The r-index /
  // sr-index variants only matter for the query binary's factory cache; they
  // aren't constructed here.)
  benchmark::RegisterBenchmark("DocListIdxBrute", BM_ConstructBruteIdx<dret::DocListIdxBrute<>>, config, data_path);

  // RMQ family — typed aliases come from bench::factories::rmq, parameterised
  // on dret::GenericStorage (each construct() takes a fresh storage).
  namespace rmq_fac = bench::factories::rmq;
  using GS = dret::GenericStorage;
  if (HasVariant(rmq_get_doc_variants, GetDocEnum::DA)) {
    using SadaIdx = rmq_fac::Idx<GS, rmq_fac::SadaCore<GS>>;
    using IlcpIdx = rmq_fac::Idx<GS, rmq_fac::IlcpCore<GS>>;
    using CilcpIdx = rmq_fac::Idx<GS, rmq_fac::CilcpCore<GS>>;
    benchmark::RegisterBenchmark("DocListSADA-DA", BM_ConstructBruteIdx<SadaIdx>, config, data_path);
    benchmark::RegisterBenchmark("DocListILCP-DA", BM_ConstructBruteIdx<IlcpIdx>, config, data_path);
    benchmark::RegisterBenchmark("DocListCILCP-DA", BM_ConstructBruteIdx<CilcpIdx>, config, data_path);
  }

  // Phase C: non-sampled SLP index, plus the three RMQ listing cores backed by
  // the bare grammar::SLP<> cache (kSLPNS). All four are parameter-free (no bs/sf
  // axis). Each gets fanned out across the bare-SLP container axis (Default =
  // sdsl::int_vector<>, DV = sdsl::dac_vector<>, VV = sdsl::vlc_vector<>); type-
  // hashes on grammar::SLP<TVars,TLens> distinguish the on-disk cache files.
  auto register_slp_ns_for_tslp = [&]<typename TSLP>(const char* suffix) {
    using IdxSLP    = bench::factories::slp_ns::Idx<GS, TSLP>;
    using GetDoc    = rmq_fac::GetDocSLP_NS<GS, TSLP>;
    using SadaCore  = rmq_fac::SadaCore<GS, GetDoc>;
    using IlcpCore  = rmq_fac::IlcpCore<GS, GetDoc>;
    using CilcpCore = rmq_fac::CilcpCore<GS, GetDoc>;
    using SadaIdx   = rmq_fac::Idx<GS, SadaCore>;
    using IlcpIdx   = rmq_fac::Idx<GS, IlcpCore>;
    using CilcpIdx  = rmq_fac::Idx<GS, CilcpCore>;

    benchmark::RegisterBenchmark(
        std::string("DocListSLP-NS") + suffix, BM_ConstructDocListIdxSLP<IdxSLP>, config);
    if (HasVariant(rmq_get_doc_variants, GetDocEnum::SLP_NS)) {
      benchmark::RegisterBenchmark(std::string("DocListSADA-SLP-NS") + suffix,
                                   BM_ConstructDocListIdxRMQCompressedNS<SadaIdx, SadaCore>, config);
      benchmark::RegisterBenchmark(std::string("DocListILCP-SLP-NS") + suffix,
                                   BM_ConstructDocListIdxRMQCompressedNS<IlcpIdx, IlcpCore>, config);
      benchmark::RegisterBenchmark(std::string("DocListCILCP-SLP-NS") + suffix,
                                   BM_ConstructDocListIdxRMQCompressedNS<CilcpIdx, CilcpCore>, config);
    }
  };

  if (HasVariant(bare_slp_variants, BareSLPVariant::IV))
    register_slp_ns_for_tslp.template operator()<bench::factories::slp_ns::BareSLP_IV>("");
  if (HasVariant(bare_slp_variants, BareSLPVariant::Raw))
    register_slp_ns_for_tslp.template operator()<bench::factories::slp_ns::BareSLP_Raw>("-Raw");
  if (HasVariant(bare_slp_variants, BareSLPVariant::DV))
    register_slp_ns_for_tslp.template operator()<bench::factories::slp_ns::BareSLP_DV>("-DV");
  if (HasVariant(bare_slp_variants, BareSLPVariant::VV))
    register_slp_ns_for_tslp.template operator()<bench::factories::slp_ns::BareSLP_VV>("-VV");

  auto block_sizes = powersOfTwo(FLAGS_min_block_size, FLAGS_max_block_size);
  auto storing_factors = powersOfTwo(FLAGS_min_storing_factor, FLAGS_max_storing_factor);
  if (!block_sizes.empty() && !storing_factors.empty()) {
    // GCDA — typed-index aliases come from bench::factories::gcda, parameterised
    // on dret::GenericStorage. Type-hashes on disk match what the factory path
    // produces because TStorage / TAlphabet / TCountIdx are identical.
    namespace gcda_fac = bench::factories::gcda;
    auto register_gcda = [&]<typename TIndex>(const char* suffix) {
      const std::string name = std::string("DocListGCDA") + suffix;
      benchmark::RegisterBenchmark(name, BM_ConstructDocListIdxGCDA<TIndex>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    };
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::Light))
      register_gcda.template operator()<gcda_fac::Idx<dret::GenericStorage>>("");
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactBP))
      register_gcda.template operator()<gcda_fac::Idx<dret::GenericStorage, gcda_fac::SLP_CompactBP>>("-CompactBP");
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactLOUDS))
      register_gcda.template operator()<gcda_fac::Idx<dret::GenericStorage, gcda_fac::SLP_CompactLOUDS>>("-CompactLOUDS");
    if (HasVariant(gcda_slp_variants, GCDASLPVariant::Combined))
      register_gcda.template operator()<gcda_fac::Idx<dret::GenericStorage, gcda_fac::SLP_Combined>>("-Combined");

    // RMQ-SLP — one register per (gcda_slp variant) requested. Each variant's
    // GetDocSLP TSLP shares its cache file with the matching DocListGCDA-*
    // build via grammar:: type-hashing on disk.
    auto register_rmq_slp_for_tslp = [&]<typename TSLP>(const char* suffix) {
      using GetDoc    = rmq_fac::GetDocSLP<GS, TSLP>;
      using SadaCore  = rmq_fac::SadaCore<GS, GetDoc>;
      using IlcpCore  = rmq_fac::IlcpCore<GS, GetDoc>;
      using CilcpCore = rmq_fac::CilcpCore<GS, GetDoc>;
      using SadaIdx   = rmq_fac::Idx<GS, SadaCore>;
      using IlcpIdx   = rmq_fac::Idx<GS, IlcpCore>;
      using CilcpIdx  = rmq_fac::Idx<GS, CilcpCore>;
      const std::string sep = (suffix[0] == '\0') ? std::string{} : std::string("-") + suffix;
      benchmark::RegisterBenchmark(std::string("DocListSADA-SLP") + sep,
                                   BM_ConstructDocListIdxRMQCompressed<SadaIdx, SadaCore>, config)
          ->ArgsProduct({block_sizes, storing_factors});
      benchmark::RegisterBenchmark(std::string("DocListILCP-SLP") + sep,
                                   BM_ConstructDocListIdxRMQCompressed<IlcpIdx, IlcpCore>, config)
          ->ArgsProduct({block_sizes, storing_factors});
      benchmark::RegisterBenchmark(std::string("DocListCILCP-SLP") + sep,
                                   BM_ConstructDocListIdxRMQCompressed<CilcpIdx, CilcpCore>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    };
    if (HasVariant(rmq_get_doc_variants, GetDocEnum::SLP)) {
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::Light))
        register_rmq_slp_for_tslp.template operator()<rmq_fac::SLP_Light>("");
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactBP))
        register_rmq_slp_for_tslp.template operator()<rmq_fac::SLP_CompactBP>("CompactBP");
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::CompactLOUDS))
        register_rmq_slp_for_tslp.template operator()<rmq_fac::SLP_CompactLOUDS>("CompactLOUDS");
      if (HasVariant(gcda_slp_variants, GCDASLPVariant::Combined))
        register_rmq_slp_for_tslp.template operator()<rmq_fac::SLP_Combined>("Combined");
    }

    // RMQ-DSLP.
    if (HasVariant(rmq_get_doc_variants, GetDocEnum::DSLP)) {
      using GetDoc    = rmq_fac::GetDocDSLP<GS>;
      using SadaCore  = rmq_fac::SadaCore<GS, GetDoc>;
      using IlcpCore  = rmq_fac::IlcpCore<GS, GetDoc>;
      using CilcpCore = rmq_fac::CilcpCore<GS, GetDoc>;
      using SadaIdx   = rmq_fac::Idx<GS, SadaCore>;
      using IlcpIdx   = rmq_fac::Idx<GS, IlcpCore>;
      using CilcpIdx  = rmq_fac::Idx<GS, CilcpCore>;
      benchmark::RegisterBenchmark("DocListSADA-DSLP", BM_ConstructDocListIdxRMQCompressed<SadaIdx, SadaCore>, config)
          ->ArgsProduct({block_sizes, storing_factors});
      benchmark::RegisterBenchmark("DocListILCP-DSLP", BM_ConstructDocListIdxRMQCompressed<IlcpIdx, IlcpCore>, config)
          ->ArgsProduct({block_sizes, storing_factors});
      benchmark::RegisterBenchmark("DocListCILCP-DSLP", BM_ConstructDocListIdxRMQCompressed<CilcpIdx, CilcpCore>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    }

    // DGCDA — one construct benchmark per requested TSLP variant.
    namespace dgcda_fac = bench::factories::dgcda;
    auto register_dgcda = [&]<typename TIndex>(const char* suffix) {
      const std::string name = std::string("DocListDGCDA") + suffix;
      benchmark::RegisterBenchmark(name, BM_ConstructDocListIdxGCDA<TIndex>, config)
          ->ArgsProduct({block_sizes, storing_factors});
    };
    if (HasVariant(dgcda_slp_variants, DGCDASLPVariant::Default))
      register_dgcda.template operator()<dgcda_fac::Idx<dret::GenericStorage>>("");
    if (HasVariant(dgcda_slp_variants, DGCDASLPVariant::OTF))
      register_dgcda.template operator()<dgcda_fac::Idx<dret::GenericStorage, dgcda_fac::SLP_OTF>>("-OTF");
    if (HasVariant(dgcda_slp_variants, DGCDASLPVariant::CRL))
      register_dgcda.template operator()<dgcda_fac::Idx<dret::GenericStorage, dgcda_fac::SLP_CRL>>("-CRL");
    if (HasVariant(dgcda_slp_variants, DGCDASLPVariant::EV))
      register_dgcda.template operator()<dgcda_fac::Idx<dret::GenericStorage, dgcda_fac::SLP_EV>>("-EV");
    if (HasVariant(dgcda_slp_variants, DGCDASLPVariant::DV))
      register_dgcda.template operator()<dgcda_fac::Idx<dret::GenericStorage, dgcda_fac::SLP_DV>>("-DV");
    if (HasVariant(dgcda_slp_variants, DGCDASLPVariant::VV))
      register_dgcda.template operator()<dgcda_fac::Idx<dret::GenericStorage, dgcda_fac::SLP_VV>>("-VV");

    // PDL build benchmarks. Typed-index aliases come from bench::factories::pdl;
    // the (codec, get-doc) pair is a compile-time axis (9 template instantiations)
    // and storage_policy is a runtime arg threaded into the constructor. Skipped
    // entirely when --pdl_variants is empty.
    namespace pdl_fac = bench::factories::pdl;
    auto register_pdl = [&]<typename TIndex>(const std::string& base) {
      for (const auto sp : pdl_storage_policies) {
        std::string name = "DocListPDL-" + base + "-" + bench::axes::EnumTraits<PDLStoragePolicy>::Name(sp);
        benchmark::RegisterBenchmark(
            name, BM_ConstructDocListIdxPDL<TIndex>, config, bench::axes::toPDLStoragePolicy(sp))
            ->ArgsProduct({block_sizes, storing_factors});
      }
    };

    for (const auto pdl_v : pdl_variants) {
      for (const auto pdl_gd : pdl_get_doc_variants) {
        std::string pair = std::string(bench::axes::EnumTraits<PDLVariant>::Name(pdl_v)) + "-" +
                           std::string(bench::axes::EnumTraits<GetDocEnum>::Name(pdl_gd));
        switch (pdl_v) {
          case PDLVariant::Plain:
            switch (pdl_gd) {
              case GetDocEnum::DA:   register_pdl.template operator()<pdl_fac::IdxPlain<GS>>(pair); break;
              case GetDocEnum::SLP:  register_pdl.template operator()<pdl_fac::IdxPlain<GS, pdl_fac::GetDocsSLP<GS>>>(pair); break;
              case GetDocEnum::DSLP: register_pdl.template operator()<pdl_fac::IdxPlain<GS, pdl_fac::GetDocsDSLP<GS>>>(pair); break;
              case GetDocEnum::SLP_NS: break;  // not a valid PDL backing
            }
            break;
          case PDLVariant::RP:
            switch (pdl_gd) {
              case GetDocEnum::DA:   register_pdl.template operator()<pdl_fac::IdxRP<GS>>(pair); break;
              case GetDocEnum::SLP:  register_pdl.template operator()<pdl_fac::IdxRP<GS, pdl_fac::GetDocsSLP<GS>>>(pair); break;
              case GetDocEnum::DSLP: register_pdl.template operator()<pdl_fac::IdxRP<GS, pdl_fac::GetDocsDSLP<GS>>>(pair); break;
              case GetDocEnum::SLP_NS: break;
            }
            break;
          case PDLVariant::BC:
            switch (pdl_gd) {
              case GetDocEnum::DA:   register_pdl.template operator()<pdl_fac::IdxBC<GS>>(pair); break;
              case GetDocEnum::SLP:  register_pdl.template operator()<pdl_fac::IdxBC<GS, pdl_fac::GetDocsSLP<GS>>>(pair); break;
              case GetDocEnum::DSLP: register_pdl.template operator()<pdl_fac::IdxBC<GS, pdl_fac::GetDocsDSLP<GS>>>(pair); break;
              case GetDocEnum::SLP_NS: break;
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
