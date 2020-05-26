//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/23/20.
//

#include <iostream>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include "factory.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_string(data_name, "data", "Data file basename.");
DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State &state) {
  for (auto _ : state) {
    std::string empty_string;
  }
}
BENCHMARK(BM_WarmUp);

// Benchmark Queries on Document Frequency Index
auto BM_QueryDocFreqIndex = [](benchmark::State &_state, const auto *_idx, const auto &_patterns) {
  for (auto _ : _state) {
    for (const auto &pattern: _patterns) {
      auto freqs = _idx->Search(pattern);
    }
  }
};

auto BM_PrintQueryDocFreqIndex =
    [](benchmark::State &_state, const auto &_idx_name, const auto *_idx, const auto &_patterns) {
      for (auto _ : _state) {
        std::ofstream out(std::string("result-") + _idx_name + ".txt");
        for (const auto &pattern: _patterns) {
          out << pattern << std::endl;
          auto freqs = _idx->Search(pattern);

          std::map<std::size_t, std::size_t> ordered_freqs{freqs.begin(), freqs.end()};
          for (const auto &item  : ordered_freqs) {
            out << "  " << item.first << ":" << item.second << std::endl;
          }
        }
      }
    };

int main(int argc, char *argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_patterns.empty() || FLAGS_data_name.empty() || FLAGS_data_dir.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  // Query patterns
  std::vector<std::string> patterns;
  {
    std::ifstream pattern_file(FLAGS_patterns.c_str(), std::ios_base::binary);
    if (!pattern_file) {
      std::cerr << "ERROR: Failed to open patterns file!" << std::endl;
      return 3;
    }

    std::string buf;
    while (std::getline(pattern_file, buf)) {
      if (buf.empty())
        continue;

      patterns.emplace_back(buf);
    }
    pattern_file.close();
  }

  sdsl::cache_config config(true, FLAGS_data_dir, FLAGS_data_name);

  Factory factory(config);

  std::vector<std::pair<const char *, Factory::Config>> index_configs = {
      {"Brute-R-Index", Factory::Config{Factory::IndexEnum::Brute}},
      {"Sada", Factory::Config{Factory::IndexEnum::SADA}},
      {"ILCP", Factory::Config{Factory::IndexEnum::ILCP}},
      {"CILCP", Factory::Config{Factory::IndexEnum::CILCP}}
  };

  std::string print_bm_prefix = "Print-";
  for (const auto &idx_config: index_configs) {
    auto index = factory.Build(idx_config.second);
    benchmark::RegisterBenchmark(idx_config.first, BM_QueryDocFreqIndex, index, patterns);
    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.first;
      benchmark::RegisterBenchmark(print_bm_name.c_str(), BM_PrintQueryDocFreqIndex, idx_config.first, index, patterns);
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
