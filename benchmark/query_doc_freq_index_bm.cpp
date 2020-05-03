//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/25/20.
//

#include <iostream>

#include <gflags/gflags.h>
#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include <rindex/r_index.hpp>

#include <doc_freq_index_brute.h>
#include <doc_freq_index_sada.h>

DEFINE_string(data, "", "Data file basename. (MANDATORY)");
DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");

const char *KEY_R_INDEX = "ri";
const char *KEY_DOC_END = "doc_end";
const char *KEY_DOC_ISA = "doc_isa";
const char *KEY_SADA_RMINQ = "sada_rminq";
const char *KEY_SADA_RMAXQ = "sada_rmaxq";

using RangeMinQuery = sdsl::rmq_succinct_sct<true>;
using RangeMaxQuery = sdsl::rmq_succinct_sct<false>;

auto BM_Doc_Freq = [](benchmark::State &_state, const auto &_idx, const auto &_patterns) {
  for (auto _ : _state) {
    for (const auto &pattern: _patterns) {
//      out << pattern << std::endl;
      auto freqs = _idx.ListWithFreq(pattern);

//      std::map<std::size_t, std::size_t> ordered_freqs{freqs.begin(), freqs.end()};
//      for (const auto &item  : ordered_freqs) {
//        out << "  " << item.first << ":" << item.second << std::endl;
//      }
    }
  }
};

//auto BM_CSA_Pos = [](benchmark::State &_state, const auto &_csa, const auto &_positions) {
//  for (auto _ : _state) {
//    for (const auto &pos : _positions) {
//      auto suffix = _csa[pos];
//    }
//  }
//};

template<typename RIndex>
class RIndexWrapper {
 public:
  RIndexWrapper(const RIndex &_r_index) : r_index_{_r_index} {}

  template<typename Pattern>
  auto Locate(const Pattern &_pattern) const {
    return const_cast<RIndex &>(r_index_).locate_all(const_cast<std::string &>(_pattern));
  }

 private :
  const RIndex &r_index_;
};

template<typename RIndex, typename CSA>
class CSAWrapper {
 public:
  CSAWrapper(const RIndex &_r_index, const CSA &_csa) : r_index_{_r_index}, csa_{_csa} {}

  template<typename Pattern>
  auto Locate(const Pattern &_pattern) const {
    return sdsl::locate(csa_, _pattern.begin(), _pattern.end());
  }

  template<typename Pattern>
  auto Search(const Pattern &_pattern) const {
    return const_cast<RIndex &>(r_index_).count(const_cast<std::string &>(_pattern));
  }

  auto operator[](std::size_t _i) const {
    return csa_[_i];
  }

 private :
  const RIndex &r_index_;
  const CSA &csa_;
};

template<typename RIndex, typename CSA>
auto MakeCSAWrapper(const RIndex &_r_index, const CSA &_csa) {
  return CSAWrapper<RIndex, CSA>(_r_index, _csa);
}

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

  //*************
  // Brute r-index
  //*************

  ri::r_index<> r_idx;
  if (!sdsl::load_from_cache(r_idx, KEY_R_INDEX, config)) {
    std::cerr << "ERROR: Failed load of R-Index!" << std::endl;
    exit(1);
  }
  RIndexWrapper<ri::r_index<>> r_idx_wrapper{r_idx};

  sdsl::bit_vector doc_endings;
  if (!sdsl::load_from_cache(doc_endings, KEY_DOC_END, config)) {
    std::cerr << "ERROR: Failed load of Documents Endings!" << std::endl;
    exit(1);
  }

  using BitVector = sdsl::sd_vector<>;

  BitVector doc_endings_compact(doc_endings);
  auto doc_endings_rank = BitVector::rank_1_type(&doc_endings_compact);
  auto doc_endings_select = BitVector::select_1_type(&doc_endings_compact);
  std::size_t doc_cnt = doc_endings_rank(doc_endings.size());

  auto idx_brute = dret::MakeDocFreqIndexBrute(r_idx_wrapper, doc_endings_rank);

  benchmark::RegisterBenchmark("Brute-RIndex", BM_Doc_Freq, idx_brute, patterns);


  //*************
  // SADA
  //*************

  sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<63>>, 30, 1000000, sdsl::text_order_sa_sampling<> > csa;
  if (!sdsl::load_from_cache(csa, conf::KEY_CSA, config)) {
    std::cerr << "ERROR: Failed load of CSA!" << std::endl;
    exit(1);
  }

  auto csa_wrapper = MakeCSAWrapper(r_idx, csa);

  vector<int_vector<> > doc_isa;
  if (!sdsl::load_from_cache(doc_isa, KEY_DOC_ISA, config)) {
    std::cerr << "ERROR: Failed load of Documents ISA!" << std::endl;
    exit(1);
  }

  RangeMinQuery range_min_query;
  if (!sdsl::load_from_cache(range_min_query, KEY_SADA_RMINQ, config)) {
    std::cerr << "ERROR: Failed load of Sada RMinQ!" << std::endl;
    exit(1);
  }

  RangeMaxQuery range_max_query;
  if (!sdsl::load_from_cache(range_max_query, KEY_SADA_RMAXQ, config)) {
    std::cerr << "ERROR: Failed load of Sada RMaxQ!" << std::endl;
    exit(1);
  }

  auto idx_sada = dret::MakeDocFreqIndexSada(csa_wrapper,
                                             doc_cnt,
                                             doc_isa,
                                             doc_endings_compact,
                                             doc_endings_rank,
                                             doc_endings_select,
                                             range_min_query,
                                             range_max_query);

  benchmark::RegisterBenchmark("Sada", BM_Doc_Freq, idx_sada, patterns);


  //*************
  // Brute csa
  //*************
//  auto idx_brute_csa = dret::MakeDocFreqIndexBrute(csa_wrapper, doc_endings_rank);
//
//  benchmark::RegisterBenchmark("Brute-CSA", BM_Doc_Freq, idx_brute_csa, patterns);

//  std::ifstream in("pos.txt");
//  std::size_t pos;
//  std::vector<std::size_t> positions;
//  while (in >> pos) {
//    positions.emplace_back(pos);
//  }
//
//  benchmark::RegisterBenchmark("CSA[]", BM_CSA_Pos, csa, positions);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}