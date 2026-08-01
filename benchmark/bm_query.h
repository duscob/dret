//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/8/26.
//

#pragma once

#include <algorithm>
#include <fstream>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include "base64.h"
#include "bm_base.h"
#include "bm_size_counters.h"

#include "dret/size_report.h"

#ifdef DRET_DOC_LIST_PROFILE
#include "dret/doc_list/search_profile.h"
#endif

DEFINE_string(pattern_code, "PLAIN", "Codification Algorithm for pattern: PLAIN, BASE64");

//~~~~~~~


void SetupDefaultCounters(benchmark::State& t_state) {
  t_state.counters["Collection_Size(bytes)"] = 0;
  t_state.counters["Size(bytes)"] = 0;
  t_state.counters["Bits_x_Symbol"] = 0;
  t_state.counters["Patterns"] = 0;
  t_state.counters["Time_x_Pattern"] = 0;
  t_state.counters["Results"] = 0;
  t_state.counters["Time_x_Result"] = 0;
}

//~~~~~~~


// Benchmark Warm-up
static void BM_WarmUp(benchmark::State& _state) {
  for (auto _ : _state) {
    std::vector<int> empty_vector(1000000, 0);
  }

  SetupDefaultCounters(_state);
}

BENCHMARK(BM_WarmUp);

//~~~~~~~


auto UpdateCounter = [](benchmark::State& t_state, auto t_n, auto t_index_size, auto t_n_patterns, auto t_n_results) {
  SetupDefaultCounters(t_state);
  t_state.counters["Collection_Size(bytes)"] = t_n;
  t_state.counters["Size(bytes)"] = t_index_size;
  t_state.counters["Bits_x_Symbol"] = t_index_size * 8.0 / t_n;
  t_state.counters["Patterns"] = t_n_patterns;
  t_state.counters["Time_x_Pattern"] =
      benchmark::Counter(t_n_patterns, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
  t_state.counters["Results"] = t_n_results;
  t_state.counters["Time_x_Result"] =
      benchmark::Counter(t_n_results, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

//~~~~~~~


auto BM_MacroQuery = [](benchmark::State& t_state, auto t_make_index, const auto& t_patterns, auto t_n) {
  auto [query, index_size, sizes] = t_make_index(t_state);

  std::size_t total = 0;

#ifdef DRET_DOC_LIST_PROFILE
  dret::search_profile().reset();
#endif

  for (auto _ : t_state) {
    total = 0;
    for (const auto& pattern : t_patterns) {
      auto results = query(pattern.decoded);
      total += results.size();
    }
  }

  UpdateCounter(t_state, t_n, index_size, t_patterns.size(), total);
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
#ifdef DRET_DOC_LIST_PROFILE
  {
    const auto& pr = dret::search_profile();
    const double q = pr.n_queries ? static_cast<double>(pr.n_queries) : 1.0;
    t_state.counters["prof_count_ns_x_q"] = pr.ns_count / q;
    t_state.counters["prof_cover_ns_x_q"] = pr.ns_cover / q;
    t_state.counters["prof_expand_ns_x_q"] = pr.ns_expand / q;
    t_state.counters["prof_combine_ns_x_q"] = pr.ns_combine / q;
    t_state.counters["prof_nodes_x_q"] = pr.n_nodes / q;
    t_state.counters["prof_rawpos_x_q"] = pr.n_raw_positions / q;
    t_state.counters["prof_docs_x_q"] = pr.n_docs / q;
  }
#endif
  {
    auto bm_name = t_state.name();
    std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
    std::replace(idx_name.begin(), idx_name.end(), '/', '_');
    dret::writeSizesJson("sizes-query-" + idx_name + ".json", sizes);
  }
};

//~~~~~~~


auto BM_MicroQuery = [](benchmark::State& t_state, auto t_make_index, const auto& t_patterns, auto t_i, auto t_n) {
  auto [query, index_size, sizes] = t_make_index(t_state);

  const auto& pattern = t_patterns[*t_i];
  std::size_t total = 0;

  for (auto _ : t_state) {
    auto results = query(pattern.decoded);
    total = results.size();
  }

  if (++*t_i == t_patterns.size()) {
    *t_i = 0;
  }

  UpdateCounter(t_state, t_n, index_size, 1, total);
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
};

//~~~~~~~


auto BM_PrintQuery =
    [](benchmark::State& t_state, auto t_make_index, const auto& t_patterns, auto t_n, bool t_is_sampled_index) {
      auto bm_name = t_state.name();
      std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
      replace(idx_name.begin(), idx_name.end(), '/', '_');
      if (t_is_sampled_index) {
        idx_name += "-" + std::to_string(t_state.range(0));
      }
      std::string output_filename = "result-" + idx_name + ".txt";

      auto [query, index_size, sizes] = t_make_index(t_state);

      std::size_t total = 0;

      for (auto _ : t_state) {
        std::ofstream out(output_filename);
        total = 0;
        for (const auto& pattern : t_patterns) {
          out << pattern.encoded << std::endl;
          auto results = query(pattern.decoded);
          total += results.size();

          PrintResults(out, results);
        }
      }

      UpdateCounter(t_state, t_n, index_size, t_patterns.size(), total);
      appendCounters(t_state, sizes);
      t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
    };

enum KeyQueryBenchmark {
  kMacro,
  kMicro,
  kPrint,
};

//~~~~~~~


using QueryBenchmarks = std::map<KeyQueryBenchmark, benchmark::internal::Benchmark*>;

struct QueryBenchmarkConfig {
  bool report_stats = false;
  int reps = 10;
  double min_time = 0.0;
  bool print_results = false;
};

//~~~~~~~


auto RegisterQueryBenchmarks = [](const auto& t_name,
                                  auto t_make_index,
                                  const auto& t_patterns,
                                  auto t_n,
                                  const QueryBenchmarkConfig& t_bm_config,
                                  bool t_is_sampled_index) {
  QueryBenchmarks bms;

  bms[kMacro] = benchmark::RegisterBenchmark(t_name, BM_MacroQuery, t_make_index, t_patterns, t_n);

  if (t_bm_config.report_stats) {
    auto statistics_min = [](const std::vector<double>& v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
    };
    auto statistics_max = [](const std::vector<double>& v) -> double {
      return *(std::max_element(std::begin(v), std::end(v)));
    };

    bms[kMacro]
        ->Name(t_name + "/macro")
        ->Repetitions(t_bm_config.reps)
        ->ComputeStatistics("min", statistics_min)
        ->ComputeStatistics("max", statistics_max)
        ->ReportAggregatesOnly();

    bms[kMicro] = benchmark::RegisterBenchmark(
                      t_name + "/micro", BM_MicroQuery, t_make_index, t_patterns, std::make_shared<int>(0), t_n)
                      ->Repetitions(t_patterns.size())
                      ->ComputeStatistics("min", statistics_min)
                      ->ComputeStatistics("max", statistics_max)
                      ->ComputeStatistics("total",
                                          [](const std::vector<double>& v) -> double {
                                            return std::accumulate(std::begin(v), std::end(v), double(0));
                                          })
                      ->ReportAggregatesOnly();
    if (t_bm_config.min_time > 0) {
      bms[kMicro]->MinTime(t_bm_config.min_time);
    } else {
      bms[kMicro]->Iterations(t_bm_config.reps);
    }
  }

  if (t_bm_config.print_results) {
    bms[kPrint] = benchmark::RegisterBenchmark(
                      "Print/" + t_name, BM_PrintQuery, t_make_index, t_patterns, t_n, t_is_sampled_index)
                      ->Iterations(1);
  }

  return bms;
};

//~~~~~~~


template <typename TSequence>
struct Pattern {
  TSequence encoded;
  TSequence decoded;

  Pattern(TSequence t_encoded, TSequence t_decoded) : encoded(std::move(t_encoded)), decoded(std::move(t_decoded)) {}
};

//~~~~~~~


enum PatternCode {
  kPlain,
  kBase64,
};

//~~~~~~~


PatternCode toPatternCode(const std::string& t_str) {
  static const std::map<std::string, PatternCode> name_to_enum = {
      {"PLAIN", kPlain},
      {"BASE64", kBase64},
  };

  return name_to_enum.at(t_str);
}

//~~~~~~~


template <typename TSequence>
auto ReadPatterns(const std::string& t_pattern_path,
                  uint8_t t_num_bytes = 1,
                  typename TSequence::value_type t_eol = '\n') {
  std::ifstream pattern_file(t_pattern_path, std::ios_base::binary);
  if (!pattern_file) {
    std::cerr << "ERROR: Failed to open patterns file! (" << t_pattern_path << ")" << std::endl;
    exit(3);
  }

  auto pattern_code = toPatternCode(FLAGS_pattern_code);

  auto decode = [pattern_code](const auto& tt_pattern) {
    switch (pattern_code) {
      case kBase64:
        return base64_decode(tt_pattern);
      default:
        return tt_pattern;
    }
  };

  std::vector<Pattern<TSequence>> patterns;
  TSequence pattern;
  typename TSequence::value_type value;
  while (pattern_file.read(reinterpret_cast<char*>(&value), t_num_bytes)) {
    if (value != t_eol) {
      pattern.push_back(value);
    } else {
      patterns.emplace_back(pattern, decode(pattern));
      pattern.clear();
    }
  }
  if (!pattern.empty()) {
    patterns.emplace_back(pattern, decode(pattern));
  }

  pattern_file.close();

  return patterns;
}
