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
#include <doc_freq_index_rmq.h>

DEFINE_string(data, "", "Data file basename. (MANDATORY)");
DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");

const char *KEY_R_INDEX = "ri";
const char *KEY_DOC_END = "doc_end";
const char *KEY_DOC_ISA = "doc_isa";
const char *KEY_SADA_RMINQ = "sada_rminq";
const char *KEY_SADA_RMAXQ = "sada_rmaxq";
const char *KEY_ILCP_BACKWARD_RUN_HEADS = "ilcp_b_run_heads";
const char *KEY_ILCP_FORWARD_RUN_HEADS = "ilcp_f_run_heads";
const char *KEY_ILCP_BACKWARD_RMQ = "ilcp_b_rmq";
const char *KEY_ILCP_FORWARD_RMQ = "ilcp_f_rmq";
const char *KEY_CILCP_BACKWARD_RUN_HEADS = "cilcp_b_run_heads";
const char *KEY_CILCP_FORWARD_RUN_HEADS = "cilcp_f_run_heads";
const char *KEY_CILCP_BACKWARD_RMQ = "cilcp_b_rmq";
const char *KEY_CILCP_FORWARD_RMQ = "cilcp_f_rmq";

using RangeMinQuery = sdsl::rmq_succinct_sct<true>;
using RangeMaxQuery = sdsl::rmq_succinct_sct<false>;

auto BM_Doc_Freq = [](benchmark::State &_state, const auto &_idx, const auto &_patterns) {
  for (auto _ : _state) {
//    std::ofstream out("result.txt");
    for (const auto &pattern: _patterns) {
//      out << pattern << std::endl;
//      std::cout << pattern << std::endl;
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

template<typename CSA, typename DocFreqIdx>
class DocFreqIdxWrapper {
 public:
  DocFreqIdxWrapper(const CSA &_csa, const DocFreqIdx &_doc_freq_idx) : csa_{_csa}, doc_freq_idx_{_doc_freq_idx} {}

  template<typename Pattern>
  auto ListWithFreq(const Pattern &_pattern) const {
    auto[sp, ep] = csa_.Search(_pattern);

    return doc_freq_idx_.ListWithFreq(sp, ep);
  }
 private:
  const CSA &csa_;
  const DocFreqIdx &doc_freq_idx_;
};

template<typename CSA, typename DocFreqIdx>
auto MakeDocFreqIdxWrapper(const CSA &_csa, const DocFreqIdx &_doc_freq_idx) {
  return DocFreqIdxWrapper<CSA, DocFreqIdx>(_csa, _doc_freq_idx);
}

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

  auto default_transform = [](auto _sp, auto _ep, auto &_stack, dret::OccurrenceSide _side) {
    _stack.emplace(_sp, _ep);
    return std::make_pair(_sp, _ep);
  };

  auto get_suffix_doc = [&csa, &doc_endings_rank](std::size_t _i, auto... params) {
    auto suffix = csa[_i];
    return std::make_pair(suffix, doc_endings_rank(suffix));
  };

  auto report_suffix_doc =
      [](auto _i, const auto &_suffix_doc, dret::OccurrenceSide _side, auto _report, auto _mark, auto... params) {
        _report(_suffix_doc.first);
        _mark(_suffix_doc.second, _side);
      };

  std::array<sdsl::bit_vector, 2> marked_docs{sdsl::bit_vector(doc_cnt, 0), sdsl::bit_vector(doc_cnt, 0)};
  auto is_marked = [&marked_docs](auto _i, const auto &_suffix_doc, dret::OccurrenceSide _side) {
    return marked_docs[_side == dret::OccurrenceSide::LEFTMOST ? 0 : 1][_suffix_doc.second];
  };

  auto mark_doc = [&marked_docs](auto _doc, dret::OccurrenceSide _side) {
    marked_docs[_side == dret::OccurrenceSide::LEFTMOST ? 0 : 1][_doc] = 1;
  };

  auto unmark_doc = [&marked_docs](auto _doc) {
    marked_docs[0][_doc] = 0;
    marked_docs[1][_doc] = 0;
  };

  auto idx_rmq_sada = dret::MakeDocFreqIndexRMQ(doc_cnt,
                                                range_min_query,
                                                range_max_query,
                                                default_transform,
                                                get_suffix_doc,
                                                report_suffix_doc,
                                                is_marked,
                                                mark_doc,
                                                unmark_doc,
                                                doc_isa,
                                                doc_endings_compact,
                                                doc_endings_rank,
                                                doc_endings_select);
  auto idx_rmq_sada_wrapper = MakeDocFreqIdxWrapper(csa_wrapper, idx_rmq_sada);

  benchmark::RegisterBenchmark("RMQ-Sada", BM_Doc_Freq, idx_rmq_sada_wrapper, patterns);


  //*************
  // ILCP
  //*************

  RangeMinQuery ilcp_rmq_left;
  if (!sdsl::load_from_cache(ilcp_rmq_left, KEY_ILCP_BACKWARD_RMQ, config)) {
    std::cerr << "ERROR: Failed load of ILCP Left RMQ!" << std::endl;
    exit(1);
  }
  RangeMinQuery ilcp_rmq_right;
  if (!sdsl::load_from_cache(ilcp_rmq_right, KEY_ILCP_FORWARD_RMQ, config)) {
    std::cerr << "ERROR: Failed load of ILCP Right RMQ!" << std::endl;
    exit(1);
  }

  sdsl::bit_vector ilcp_run_heads_left;
  if (!sdsl::load_from_cache(ilcp_run_heads_left, KEY_ILCP_BACKWARD_RUN_HEADS, config)) {
    std::cerr << "ERROR: Failed load of ILCP Left Run Heads!" << std::endl;
    exit(1);
  }
  sdsl::bit_vector ilcp_run_heads_right;
  if (!sdsl::load_from_cache(ilcp_run_heads_right, KEY_ILCP_FORWARD_RUN_HEADS, config)) {
    std::cerr << "ERROR: Failed load of ILCP Right Run Heads!" << std::endl;
    exit(1);
  }

  BitVector ilcp_run_heads_compact[2] = {ilcp_run_heads_left, ilcp_run_heads_right};
  BitVector::rank_1_type ilcp_run_heads_compact_rank[2] =
      {BitVector::rank_1_type(&ilcp_run_heads_compact[0]), BitVector::rank_1_type(&ilcp_run_heads_compact[1])};
  BitVector::select_1_type ilcp_run_heads_compact_select[2] =
      {BitVector::select_1_type(&ilcp_run_heads_compact[0]), BitVector::select_1_type(&ilcp_run_heads_compact[1])};
  std::cout << "# Runs ILCP: " << ilcp_run_heads_compact_rank[0](csa.size()) << std::endl;

  auto ilcp_transform_range = /*[](auto &ilcp_run_heads_compact_rank) {
    return*/ [&ilcp_run_heads_compact_rank, &ilcp_run_heads_compact_select](auto _sp, auto _ep, auto &_stack, dret::OccurrenceSide _side) {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

//    std::cout << "O  " << _sp << "-" << _ep << std::endl;
    auto p = std::make_pair(ilcp_run_heads_compact_rank[direction](_sp + 1) - 1,
                          ilcp_run_heads_compact_rank[direction](_ep + 1) - 1);
    _stack.emplace(p);
//    std::cout << "C  " << p.first << "-" << p.second << std::endl;
//    std::cout << "R  " << ilcp_run_heads_compact_select[direction](p.first + 1) << "-" << ilcp_run_heads_compact_select[direction](p.second + 1) << std::endl;

    return p;
  };
//  };

  auto ilcp_get_suffix_doc = /*[&get_suffix_doc](auto &_ilcp_run_heads_compact_select) {
    return*/ [&ilcp_run_heads_compact_select, &get_suffix_doc](
      std::size_t _i, dret::OccurrenceSide _side, std::size_t _sp, auto... params) {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

//    std::cout << "  k': " << _i << std::endl;
    auto pos = std::max(_sp, ilcp_run_heads_compact_select[direction](_i + 1));
//    std::cout << "  k : " << pos << std::endl;

    auto v = get_suffix_doc(pos);
//    std::cout << "    d: " << v.second << std::endl;

    return v;
  };
//  };

  auto ilcp_report_suffix_doc = [&ilcp_run_heads_compact_select, &report_suffix_doc, &get_suffix_doc](
      auto _i,
      const auto &_suffix_doc,
      dret::OccurrenceSide _side,
      auto _report,
      auto _mark,
      auto _sp,
      auto _ep,
      auto... params) {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

    report_suffix_doc(_i, _suffix_doc, _side, _report, _mark);

    auto b = std::max(_sp, ilcp_run_heads_compact_select[direction](_i + 1)) + 1;
    auto e = std::min(_ep, ilcp_run_heads_compact_select[direction](_i + 2) - 1);

    for (auto i = b; i <= e; ++i) {
      auto suffix_doc = get_suffix_doc(i);
      report_suffix_doc(i, suffix_doc, _side, _report, _mark);
    }
  };

  auto idx_ilcp = dret::MakeDocFreqIndexRMQ(doc_cnt,
                                            ilcp_rmq_left,
                                            ilcp_rmq_right,
                                            ilcp_transform_range,//(ilcp_run_heads_compact_rank),
                                            ilcp_get_suffix_doc,//(ilcp_run_heads_compact_select),
                                            ilcp_report_suffix_doc,
                                            is_marked,
                                            mark_doc,
                                            unmark_doc,
                                            doc_isa,
                                            doc_endings_compact,
                                            doc_endings_rank,
                                            doc_endings_select);

  auto idx_ilcp_wrapper = MakeDocFreqIdxWrapper(csa_wrapper, idx_ilcp);

  benchmark::RegisterBenchmark("ILCP", BM_Doc_Freq, idx_ilcp_wrapper, patterns);


  //*************
  // CILCP
  //*************

  RangeMinQuery cilcp_rmq_left;
  if (!sdsl::load_from_cache(cilcp_rmq_left, KEY_CILCP_BACKWARD_RMQ, config)) {
    std::cerr << "ERROR: Failed load of CILCP Left RMQ!" << std::endl;
    exit(1);
  }
  RangeMinQuery cilcp_rmq_right;
  if (!sdsl::load_from_cache(cilcp_rmq_right, KEY_CILCP_FORWARD_RMQ, config)) {
    std::cerr << "ERROR: Failed load of CILCP Right RMQ!" << std::endl;
    exit(1);
  }

  sdsl::bit_vector cilcp_run_heads_left;
  if (!sdsl::load_from_cache(cilcp_run_heads_left, KEY_CILCP_BACKWARD_RUN_HEADS, config)) {
    std::cerr << "ERROR: Failed load of CILCP Left Run Heads!" << std::endl;
    exit(1);
  }
  sdsl::bit_vector cilcp_run_heads_right;
  if (!sdsl::load_from_cache(cilcp_run_heads_right, KEY_CILCP_FORWARD_RUN_HEADS, config)) {
    std::cerr << "ERROR: Failed load of CILCP Right Run Heads!" << std::endl;
    exit(1);
  }

  BitVector cilcp_run_heads_compact[2] = {cilcp_run_heads_left, cilcp_run_heads_right};
  BitVector::rank_1_type cilcp_run_heads_compact_rank[2] =
      {BitVector::rank_1_type(&cilcp_run_heads_compact[0]), BitVector::rank_1_type(&cilcp_run_heads_compact[1])};
  BitVector::select_1_type cilcp_run_heads_compact_select[2] =
      {BitVector::select_1_type(&cilcp_run_heads_compact[0]), BitVector::select_1_type(&cilcp_run_heads_compact[1])};

  auto cilcp_transform_range = /*[](auto &ilcp_run_heads_compact_rank) {
    return*/ [&cilcp_run_heads_compact_rank, &cilcp_run_heads_compact_select](auto _sp, auto _ep, auto &_stack, dret::OccurrenceSide _side) {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

//    std::cout << "O  " << _sp << " - " << _ep << std::endl;
    auto p = std::make_pair(cilcp_run_heads_compact_rank[direction](_sp + 1) - 1,
                            cilcp_run_heads_compact_rank[direction](_ep + 1) - 1);

    if (p.first < p.second) {
      if (_side == dret::OccurrenceSide::LEFTMOST) {
        _stack.emplace(p.second, p.second);
        _stack.emplace(p.first, p.second - 1);
      }
      else if (_side == dret::OccurrenceSide::RIGHTMOST) {
        _stack.emplace(p.first, p.first);
        _stack.emplace(p.first + 1, p.second);
      }
    }
    else {
      _stack.emplace(p);
    }
//    std::cout << "C  " << p.first << " - " << p.second << std::endl;
//    std::cout << "R  " << cilcp_run_heads_compact_select[direction](p.first + 1) << " - " << cilcp_run_heads_compact_select[direction](p.second + 1) << " - " << cilcp_run_heads_compact_select[direction](p.second + 2) << std::endl;

    return p;
  };
//  };

  auto cilcp_get_suffix_doc = /*[&get_suffix_doc](auto &_ilcp_run_heads_compact_select) {
    return*/ [&cilcp_run_heads_compact_select, &get_suffix_doc](
      std::size_t _i, dret::OccurrenceSide _side, std::size_t _sp, std::size_t _ep, auto... params) {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

//    std::cout << "  k': " << _i << std::endl;
    auto pos = direction == 0 ? std::max(_sp, cilcp_run_heads_compact_select[direction](_i + 1)) :
        std::min(_ep, cilcp_run_heads_compact_select[direction](_i + 2) - 1);
//    std::cout << "  k : " << pos << std::endl;

    auto v = get_suffix_doc(pos);
//    std::cout << "    d: " << v.second << std::endl;

    return v;
  };
//  };

  auto cilcp_report_suffix_doc = [&cilcp_run_heads_compact_select, &report_suffix_doc, &get_suffix_doc](
      auto _i,
      const auto &_suffix_doc,
      dret::OccurrenceSide _side,
      auto _report,
      auto _mark,
      auto _sp,
      auto _ep,
      auto... params) {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

    report_suffix_doc(_i, _suffix_doc, _side, _report, _mark);

    auto b = std::max(_sp, cilcp_run_heads_compact_select[direction](_i + 1)) + 1 - direction;
    auto e = std::min(_ep, cilcp_run_heads_compact_select[direction](_i + 2) - 1) - direction;

    for (auto i = b; i <= e; ++i) {
      auto suffix_doc = get_suffix_doc(i);
      if (_suffix_doc.second == suffix_doc.second) break;

      report_suffix_doc(i, suffix_doc, _side, _report, _mark);
    }
  };

  auto idx_cilcp = dret::MakeDocFreqIndexRMQ(doc_cnt,
                                            cilcp_rmq_left,
                                            cilcp_rmq_right,
                                            cilcp_transform_range,//(cilcp_run_heads_compact_rank),
                                            cilcp_get_suffix_doc,//(cilcp_run_heads_compact_select),
                                            cilcp_report_suffix_doc,
                                            is_marked,
                                            mark_doc,
                                            unmark_doc,
                                            doc_isa,
                                            doc_endings_compact,
                                            doc_endings_rank,
                                            doc_endings_select);

  auto idx_cilcp_wrapper = MakeDocFreqIdxWrapper(csa_wrapper, idx_cilcp);

  benchmark::RegisterBenchmark("CILCP", BM_Doc_Freq, idx_cilcp_wrapper, patterns);
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