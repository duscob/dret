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
#include "dret/index_base.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_string(sa_algo, "SDSL_SE_SAIS", "Suffix Array Algorithm: SDSL_SE_SAIS, SDSL_LIBDIVSUFSORT, BIG_BWT");
DEFINE_int32(doc_delim, 3, "Document delimiter.");

//~~~~~~~


void SetupCommonCounters(benchmark::State& t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["r"] = 0;
  t_state.counters["s"] = 0;
  t_state.counters["r'"] = 0;
  t_state.counters["mr'"] = 0;
}

//~~~~~~~


/// Build text representation to use with SDSL functionalities.
void BM_BuildText(benchmark::State& t_state, dret::Config t_config, const std::string& t_data_path) {
  std::size_t n;
  for (auto _ : t_state) {
    if (!cache_file_exists(sdsl::conf::KEY_TEXT, t_config)) {
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
    }
  }

  SetupCommonCounters(t_state);
  t_state.counters["n"] = n;
};

//~~~~~~~


/// Build Suffix Array
void BM_BuildSA(benchmark::State& t_state, dret::Config t_config) {
  for (auto _ : t_state) {
    if (!cache_file_exists(sdsl::conf::KEY_SA, t_config)) {
      auto event = sdsl::memory_monitor::event("SA");

      // Use SDSL functionality to build the SA
      sdsl::construct_config().byte_algo_sa =
          t_config.sa_algo == sri::SDSL_LIBDIVSUFSORT ? sdsl::LIBDIVSUFSORT : sdsl::SE_SAIS;
      sdsl::construct_sa<8>(t_config);
    }
  }

  SetupCommonCounters(t_state);
}

//~~~~~~~


/// Build document endings
void BM_BuildDocEndings(benchmark::State& t_state, dret::Config t_config) {
  for (auto _ : t_state) {
    if (!cache_file_exists(dret::conf::KEY_DOC_END, t_config)) {
      auto event = sdsl::memory_monitor::event("Doc Ends");
      dret::ConstructDocEnd(t_config);
    }
  }

  SetupCommonCounters(t_state);
}

//~~~~~~~


/// Build document array
void BM_BuildDocArray(benchmark::State& t_state, dret::Config t_config) {
  for (auto _ : t_state) {
    if (!cache_file_exists(dret::conf::KEY_DA, t_config)) {
      auto event = sdsl::memory_monitor::event("DA");
      dret::ConstructDocArray(t_config);
    }
  }

  SetupCommonCounters(t_state);
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

  dret::Config config(data_path, std::filesystem::current_path(), sri::toSAAlgo(FLAGS_sa_algo));

  benchmark::RegisterBenchmark("BuildText", BM_BuildText, config, data_path);
  benchmark::RegisterBenchmark("BuildSA", BM_BuildSA, config);
  benchmark::RegisterBenchmark("BuildDocEndings", BM_BuildDocEndings, config);
  benchmark::RegisterBenchmark("BuildDA", BM_BuildDocArray, config);

  benchmark::Initialize(&argc, argv);

  sdsl::memory_monitor::start();
  benchmark::RunSpecifiedBenchmarks();
  sdsl::memory_monitor::stop();
  {
    std::ofstream ofs("construction-common-items.html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-common-items.json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
    ofs.close();
  }

  return 0;
}
