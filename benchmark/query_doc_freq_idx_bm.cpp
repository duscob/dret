//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/25/20.
//

#include <iostream>

#include <gflags/gflags.h>
#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <rindex/r_index.hpp>

#include <df_brute.h>

DEFINE_string(data, "", "Data file basename. (MANDATORY)");
DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");

const char *KEY_R_INDEX = "ri";
const char *KEY_DOC_END = "doc_end";

auto BM_Doc_Freq = [](benchmark::State &_state, const auto &_idx, const auto &_patterns) {
  for (auto _ : _state) {
    for (const auto &pattern: _patterns) {
      _idx.ListWithFreq(pattern);
    }
  }
};


template<typename CSA>
class CSAWrapper {
 public:
  CSAWrapper(const CSA &_csa): csa_{_csa} {}

  template<typename Pattern>
  auto Locate(const Pattern &_pattern) const {
    return const_cast<CSA &>(csa_).locate_all(const_cast<std::string &>(_pattern));
  }

 private :
  const CSA &csa_;
};


int main(int argc, char *argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty() || FLAGS_patterns.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  // Query patterns
  std::vector<std::string> patterns;
  {
    std::ifstream pattern_file(FLAGS_patterns.c_str(), std::ios_base::binary);
    if (!pattern_file) {
      std::cerr << "ERROR: Failed opening of patterns file!" << std::endl;
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

  sdsl::cache_config config(false, ".", FLAGS_data);

  ri::r_index<> r_idx;
  if (!sdsl::load_from_cache(r_idx, KEY_R_INDEX, config)) {
    std::cerr << "ERROR: Failed load of R-Index!" << std::endl;
    exit(1);
  }
  CSAWrapper<ri::r_index<>> r_idx_wrapper{r_idx};

  sdsl::bit_vector doc_endings;
  if (!sdsl::load_from_cache(doc_endings, KEY_DOC_END, config)) {
    std::cerr << "ERROR: Failed load of Documents Endings!" << std::endl;
    exit(1);
  }

  using BitVector = sdsl::sd_vector<>;

  BitVector doc_endings_compact(doc_endings);
  auto doc_endings_rank = BitVector::rank_1_type(&doc_endings_compact);

  auto idx = dret::MakeDFIdxBrute(r_idx_wrapper, doc_endings_rank);

  benchmark::RegisterBenchmark("Brute", BM_Doc_Freq, idx, patterns);


  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}