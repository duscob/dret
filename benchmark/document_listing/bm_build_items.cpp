//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/8/21.
//

#include <algorithm>
#include <iostream>

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/vlc_vector.hpp>

#include "dret/basic_slp_span_length.h"
#include "dret/construct_base.h"
#include "dret/differential_light_slp.h"
#include "dret/doc_list_index_brute.h"
#include "dret/doc_list_sampled_tree_dgcda.h"
#include "dret/doc_list_sampled_tree_gcda.h"
#include "dret/index_base.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_int32(data_width, 8, "Data width in bits: 8, 16, 32 or 64");
DEFINE_string(sa_algo, "SDSL_SE_SAIS", "Suffix Array Algorithm: SDSL_SE_SAIS, SDSL_LIBDIVSUFSORT, BIG_BWT");
DEFINE_int32(doc_delim, 0, "Document delimiter.");

DEFINE_int32(min_block_size, 512, "Minimum block size (power of 2).");
DEFINE_int32(max_block_size, 512, "Maximum block size (power of 2).");
DEFINE_int32(min_storing_factor, 4, "Minimum storing factor (power of 2).");
DEFINE_int32(max_storing_factor, 4, "Maximum storing factor (power of 2).");

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

  dret::Config config(data_path,
                      std::filesystem::current_path(),
                      sri::toSAAlgo(FLAGS_sa_algo),
                      false,
                      FLAGS_data_width,
                      FLAGS_doc_delim);

  benchmark::RegisterBenchmark("DocListIdxBrute", BM_ConstructBruteIdx<dret::DocListIdxBrute<>>, config, data_path);

  auto block_sizes = powersOfTwo(FLAGS_min_block_size, FLAGS_max_block_size);
  auto storing_factors = powersOfTwo(FLAGS_min_storing_factor, FLAGS_max_storing_factor);
  if (!block_sizes.empty() && !storing_factors.empty()) {
    benchmark::RegisterBenchmark(
        "DocListGCDA", BM_ConstructDocListIdxGCDA<dret::gcda::DocListIdxGCDA<>>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    benchmark::RegisterBenchmark(
        "DocListDGCDA", BM_ConstructDocListIdxGCDA<dret::dgcda::DocListIdxDGCDA<>>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    using DGCDA_OTF = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
        dret::DifferentialLightSLP<dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-OTF", BM_ConstructDocListIdxGCDA<DGCDA_OTF>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    using DGCDA_CRL = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
        dret::DifferentialLightSLP<dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-CRL", BM_ConstructDocListIdxGCDA<DGCDA_CRL>, config)
        ->ArgsProduct({block_sizes, storing_factors});

    // Compressed int-vector variants: vary TRoots/TSpanSums/TSamples in lockstep;
    // TSampleRootsPos stays at the class default sdsl::enc_vector<>.
    using DGCDA_EV = dret::dgcda::DocListIdxDGCDA<
        dret::GenericStorage,
        dret::Alphabet<>,
        sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
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
        sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
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
        sri::SrIdxGeneric<sri::SrIndexValidArea<dret::GenericStorage, dret::Alphabet<>>, 16>,
        dret::DifferentialLightSLP<grammar::SLP<>,
                                    grammar::SampledSLP<>,
                                    sdsl::vlc_vector<>,
                                    sdsl::vlc_vector<>,
                                    sdsl::vlc_vector<>>>;
    benchmark::RegisterBenchmark(
        "DocListDGCDA-VV", BM_ConstructDocListIdxGCDA<DGCDA_VV>, config)
        ->ArgsProduct({block_sizes, storing_factors});
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
