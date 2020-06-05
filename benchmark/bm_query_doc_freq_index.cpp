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
static void BM_WarmUp(benchmark::State &_state) {
  for (auto _ : _state) {
    std::string empty_string;
  }

  _state.counters["Size(bytes)"] = 0;
  _state.counters["Bits_x_Symbol"] = 0;
  _state.counters["Patterns"] = 0;
  _state.counters["Time_x_Pattern"] = 0;
}
BENCHMARK(BM_WarmUp);

// Benchmark Queries on Document Frequency Index
auto BM_QueryDocFreqIndex =
    [](benchmark::State &_state, const auto &_idx, const auto &_patterns, auto _seq_size) {
      for (auto _ : _state) {
        for (const auto &pattern: _patterns) {
          auto freqs = _idx.first->Search(pattern);
        }
      }

      _state.counters["Size(bytes)"] = _idx.second;
      _state.counters["Bits_x_Symbol"] = _idx.second * 8.0 / _seq_size;
      _state.counters["Patterns"] = _patterns.size();
      _state.counters["Time_x_Pattern"] = benchmark::Counter(
          _patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    };

auto BM_PrintQueryDocFreqIndex =
    [](benchmark::State &_state, const auto &_idx_name, const auto &_idx, const auto &_patterns, auto _seq_size) {
      for (auto _ : _state) {
        std::ofstream out(std::string("result-") + _idx_name + ".txt");
        for (const auto &pattern: _patterns) {
          out << pattern << std::endl;
          auto freqs = _idx.first->Search(pattern);

          std::map<std::size_t, std::size_t> ordered_freqs{freqs.begin(), freqs.end()};
          for (const auto &item  : ordered_freqs) {
            out << "  " << item.first << ":" << item.second << std::endl;
          }
        }
      }

      _state.counters["Size(bytes)"] = _idx.second;
      _state.counters["Bits_x_Symbol"] = _idx.second * 8.0 / _seq_size;
      _state.counters["Patterns"] = _patterns.size();
      _state.counters["Time_x_Pattern"] = benchmark::Counter(
          _patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
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
      {"WT", Factory::Config{Factory::IndexEnum::WT}},
      {"Sada", Factory::Config{Factory::IndexEnum::SADA_ISAs}},
      {"Sada-WT-DA", Factory::Config{Factory::IndexEnum::SADA_WT_DA}},
      {"ILCP", Factory::Config{Factory::IndexEnum::ILCP_ISAs}},
      {"ILCP-WT-DA", Factory::Config{Factory::IndexEnum::ILCP_WT_DA}},
      {"CILCP", Factory::Config{Factory::IndexEnum::CILCP_ISAs}},
      {"CILCP-WT-DA", Factory::Config{Factory::IndexEnum::CILCP_WT_DA}},
      {"GCDA", Factory::Config{Factory::IndexEnum::GCDA_ISAs}},
      {"GCDA-WT-DA", Factory::Config{Factory::IndexEnum::GCDA_WT_DA}},
      {"Full-GCDA", Factory::Config{Factory::IndexEnum::FULL_GCDA}}
  };

  std::string print_bm_prefix = "Print-";
  for (const auto &idx_config: index_configs) {
    auto index = factory.Build(idx_config.second);
    benchmark::RegisterBenchmark(idx_config.first, BM_QueryDocFreqIndex, index, patterns, factory.SequenceSize());
    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.first;
      benchmark::RegisterBenchmark(
          print_bm_name.c_str(), BM_PrintQueryDocFreqIndex, idx_config.first, index, patterns, factory.SequenceSize());
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
