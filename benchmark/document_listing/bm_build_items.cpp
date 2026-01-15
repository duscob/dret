//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/8/21.
//

#include <algorithm>
#include <iostream>

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>

#include "dret/construct_base.h"
#include "dret/doc_list_index_brute.h"
#include "dret/doc_list_sampled_tree.h"
#include "dret/index_base.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_int32(data_width, 8, "Data width in bits: 8, 16, 32 or 64");
DEFINE_string(sa_algo, "SDSL_SE_SAIS", "Suffix Array Algorithm: SDSL_SE_SAIS, SDSL_LIBDIVSUFSORT, BIG_BWT");
DEFINE_int32(doc_delim, 0, "Document delimiter.");

//~~~~~~~


void SetupCommonCounters(benchmark::State& t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["bs"] = 0;  // Block size
  t_state.counters["sf"] = 0;  // Storing factor
}

//~~~~~~~


/// Build text representation to use with SDSL functionalities.
void BM_BuildText(benchmark::State& t_state, dret::Config t_config, const std::string& t_data_path) {
  std::size_t n = 0;
  for (auto _ : t_state) {
    if (cache_file_exists(sdsl::conf::KEY_TEXT, t_config)) {
      continue;
    }

    auto event = sdsl::memory_monitor::event("Text");
    sdsl::int_vector<8> text;
    {
      // Load input text by streaming from disc
      std::string input;
      {
        std::ifstream fs(t_data_path);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
      }

      // Construct text representation for SDSL use.
      n = input.size();
      text.resize(input.size() + 1);

      std::replace_copy(input.begin(), input.end(), text.begin(), '\0', static_cast<char>(FLAGS_doc_delim));

      text[text.size() - 1] = 0;  // Append symbol zero at the end
    }

    sdsl::store_to_cache(text, sdsl::conf::KEY_TEXT, t_config);
    //    sdsl::util::clear(text);

    t_state.counters["n"] = n;
  }

  SetupCommonCounters(t_state);
};

//~~~~~~~


/// Build Suffix Array
void BM_BuildSA(benchmark::State& t_state, dret::Config t_config) {
  for (auto _ : t_state) {
    if (cache_file_exists(sdsl::conf::KEY_SA, t_config)) {
      continue;
    }

    auto event = sdsl::memory_monitor::event("SA");

    // Use SDSL functionality to build the SA
    sdsl::construct_config().byte_algo_sa =
        t_config.sa_algo == sri::SDSL_LIBDIVSUFSORT ? sdsl::LIBDIVSUFSORT : sdsl::SE_SAIS;
    sdsl::construct_sa<8>(t_config);
  }

  SetupCommonCounters(t_state);
}

//~~~~~~~


/// Build document endings
void BM_BuildDocEndings(benchmark::State& t_state, dret::Config t_config) {
  for (auto _ : t_state) {
    if (cache_file_exists(dret::conf::KEY_DOC_END, t_config)) {
      continue;
    }

    auto event = sdsl::memory_monitor::event("Doc Ends");
    dret::ConstructDocEnd(t_config);
  }

  SetupCommonCounters(t_state);
}

//~~~~~~~


/// Build document array
void BM_BuildDocArray(benchmark::State& t_state, dret::Config t_config) {
  for (auto _ : t_state) {
    if (cache_file_exists(dret::conf::KEY_DA, t_config)) {
      continue;
    }

    auto event = sdsl::memory_monitor::event("DA");
    dret::ConstructDocArray(t_config);
  }

  SetupCommonCounters(t_state);
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
    ofs.close();
  }
  {
    std::ofstream ofs("construction-" + idx_name + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
    ofs.close();
  }

  SetupCommonCounters(t_state);
  {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
    t_state.counters["n"] = buf.size();
  }
  // {
  //   sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(sri::conf::KEY_BWT_RUN_FIRST, t_config));
  //   t_state.counters["r"] = buf.size();
  // }
}

//~~~~~~~


void BM_ConstructGCDA(benchmark::State& t_state, sri::Config t_config, const std::string& t_data_path) {
  std::size_t block_size = t_state.range(0);
  std::size_t storing_factor = t_state.range(1);

  dret::GCDA<> index(block_size, storing_factor);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    dret::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  {
    std::ofstream ofs("construction-" + idx_name + "-" + std::to_string(block_size) + "-"
                      + std::to_string(storing_factor) + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-" + idx_name + "-" + std::to_string(block_size) + "-"
                      + std::to_string(storing_factor) + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
    ofs.close();
  }

  SetupCommonCounters(t_state);
  t_state.counters["bs"] = block_size;
  t_state.counters["sf"] = storing_factor;
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

  // benchmark::RegisterBenchmark("BuildText", BM_BuildText, config, data_path);
  // benchmark::RegisterBenchmark("BuildSA", BM_BuildSA, config);
  // benchmark::RegisterBenchmark("BuildDocEndings", BM_BuildDocEndings, config);
  // benchmark::RegisterBenchmark("BuildDA", BM_BuildDocArray, config);

  benchmark::RegisterBenchmark("SrIndex-Brute", BM_ConstructBruteIdx<dret::DocListIdxBrute<>>, config, data_path);

  benchmark::RegisterBenchmark("GCDA", BM_ConstructGCDA, config, data_path)->ArgsProduct({{512}, {4}});

  benchmark::Initialize(&argc, argv);

  // sdsl::memory_monitor::start();
  benchmark::RunSpecifiedBenchmarks();
  // sdsl::memory_monitor::stop();
  // {
  //   std::ofstream ofs("construction-common-items.html");
  //   sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  //   ofs.close();
  // }
  // {
  //   std::ofstream ofs("construction-common-items.json");
  //   sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  //   ofs.close();
  // }

  return 0;
}
