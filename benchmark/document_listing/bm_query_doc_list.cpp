//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/3/21.
//

#include <bitset>
#include <functional>
#include <iostream>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>

#include "factory.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_string(data_name, "data", "Data file basename.");
DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

static void SetupDefaultCounters(benchmark::State& t_state) {
  t_state.counters["Size(bytes)"] = 0;
  t_state.counters["Bits_x_Symbol"] = 0;
  t_state.counters["Patterns"] = 0;
  t_state.counters["Time_x_Pattern"] = 0;
}

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State& t_state) {
  for (auto _ : t_state) {
    std::vector<int> empty_vector(1000000, 0);
  }

  SetupDefaultCounters(t_state);
}

BENCHMARK(BM_WarmUp);

class DocListResult {
 public:
  virtual ~DocListResult() = default;

  virtual void operator()(std::size_t t_doc) = 0;

  virtual void operator()() {}

  virtual void Print(std::ostream& t_os) const = 0;
};

class DocListResultVector : public DocListResult {
 public:
  void operator()(std::size_t t_doc) override {
    result_.emplace_back(t_doc);
  }

  void operator()() override {
    sort(result_.begin(), result_.end());
    result_.erase(unique(result_.begin(), result_.end()), result_.end());
  }

  void Print(std::ostream& t_os) const override {
    for (const auto& item : result_) {
      t_os << item << '\n';
    }
  }

 private:
  std::vector<std::size_t> result_;
};

class DocListResultBitvector : public DocListResult {
 public:
  explicit DocListResultBitvector(std::size_t t_size) : result_(t_size) {}

  void operator()(std::size_t t_doc) override {
    result_[t_doc] = true;
  }

  void Print(std::ostream& t_os) const override {
    for (int i = 0; i < result_.size(); ++i) {
      if (result_[i])
        t_os << i << '\n';
    }
  }

 private:
  sdsl::bit_vector result_;
};

class DocListResultVectorBool : public DocListResult {
 public:
  explicit DocListResultVectorBool(std::size_t t_size) : result_(t_size) {}

  void operator()(std::size_t t_doc) override {
    result_[t_doc] = true;
  }

  void Print(std::ostream& t_os) const override {
    for (int i = 0; i < result_.size(); ++i) {
      if (result_[i])
        t_os << i << '\n';
    }
  }

 private:
  std::vector<bool> result_;
};

// Benchmark Queries on Document Listing Index
auto BM_QueryDocList = [](benchmark::State& t_state,  //
                          auto* t_factory,
                          const auto& t_config,
                          const auto& t_patterns,
                          auto t_seq_size) {
  auto [idx, idx_size] = t_factory->Make(t_config.config);

  for (auto _ : t_state) {
    for (const auto& pattern : t_patterns) {
      auto result = t_config.create_result();
      idx->Search(pattern, std::ref(*result));

      (*result)();

      delete result;
    }
  }

  SetupDefaultCounters(t_state);
  t_state.counters["Size(bytes)"] = idx_size;
  t_state.counters["Bits_x_Symbol"] = idx_size * 8.0 / t_seq_size;
  t_state.counters["Patterns"] = t_patterns.size();
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

auto BM_PrintQueryDocList = [](benchmark::State& t_state,
                               const auto& t_idx_name,
                               auto* t_factory,
                               const auto& t_config,
                               const auto& t_patterns) {
  auto [idx, idx_size] = t_factory->Make(t_config.config);

  for (auto _ : t_state) {
    std::ofstream out(std::string("result-") + t_idx_name + ".txt");
    for (const auto& pattern : t_patterns) {
      out << pattern << std::endl;

      auto result = t_config.create_result();
      idx->Search(pattern, std::ref(*result));

      (*result)();

      result->Print(out);

      delete result;
    }
  }

  SetupDefaultCounters(t_state);
};

int main(int argc, char* argv[]) {
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

  // Benchmarks configs
  sri::Config config(FLAGS_data_name, FLAGS_data_dir, sri::SDSL_LIBDIVSUFSORT, true);

  Factory<> factory(config);

  struct Config {
    const char* name;
    Factory<>::Config config;
    std::function<DocListResult*()> create_result;
  };

  auto create_result_vector = []() {
    return new DocListResultVector();
  };
  auto create_result_bitvector = [n_docs = factory.NDocs()]() {
    return new DocListResultBitvector(n_docs);
  };
  auto create_result_vector_bool = [n_docs = factory.NDocs()]() {
    return new DocListResultVectorBool(n_docs);
  };

  std::vector<Config> index_configs = {
      {"Brute-R-Index-V", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, create_result_vector},
      {"Brute-R-Index-Bv", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, create_result_bitvector},
      {"Brute-R-Index-VB", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, create_result_vector_bool},
      {"Brute-SR-Index-V", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, create_result_vector},
      {"Brute-SR-Index-Bv", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, create_result_bitvector},
      {"Brute-SR-Index-VB", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, create_result_vector_bool},
  };

  auto seq_size = factory.SequenceSize();
  std::string print_bm_prefix = "Print-";
  for (const auto& idx_config : index_configs) {
    benchmark::RegisterBenchmark(idx_config.name, BM_QueryDocList, &factory, idx_config, patterns, seq_size);

    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.name;
      benchmark::RegisterBenchmark(
          print_bm_name.c_str(), BM_PrintQueryDocList, idx_config.name, &factory, idx_config, patterns);
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  //  benchmark::Shutdown();

  return 0;
}
