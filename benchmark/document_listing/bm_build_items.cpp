//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/8/21.
//

#include <algorithm>
#include <iostream>

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>

#include "dret/construct_base.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_bool(rebuild, false, "Rebuild all the items.");
DEFINE_bool(sais, true, "SE_SAIS or LIBDIVSUFSORT algorithm for Suffix Array construction.");

const char kDocDelimiter = '\3';

void SetupCommonCounters(benchmark::State& t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["r"] = 0;
  t_state.counters["s"] = 0;
  t_state.counters["r'"] = 0;
  t_state.counters["mr'"] = 0;
}

/// Build document endings
void BM_BuildDocEndings(benchmark::State& t_state, sdsl::cache_config* t_config) {

  for (auto _ : t_state) {
    auto event = sdsl::memory_monitor::event("Doc Ends");
    dret::ConstructDocEnd(*t_config);
  }

  SetupCommonCounters(t_state);
};

int main(int argc, char** argv) {
  gflags::SetUsageMessage("This program calculates the ri items for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::construct_config().byte_algo_sa =
      FLAGS_sais ? sdsl::SE_SAIS
                 : sdsl::LIBDIVSUFSORT;  // or LIBDIVSUFSORT for less space-efficient but faster construction

  std::string data_path = FLAGS_data;

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

  if (!cache_file_exists(dret::conf::KEY_DOC_END, config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildDocEndings", BM_BuildDocEndings, &config);
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
