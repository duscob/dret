//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/4/20.
//

#include <iostream>
#include <algorithm>
#include <cstdio>

#include <boost/filesystem.hpp>

#include <gflags/gflags.h>

#include <sdsl/io.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/rmq_succinct_sada.hpp>

#include <rindex/r_index.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

//#include <dret/diff_slp.h>
#include <dret/doc_freq_index_sada.h>

#include "definitions.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_bool(sais, true, "SE_SAIS or LIBDIVSUFSORT algorithm for Suffix Array construction.");

using namespace sdsl;
using namespace std;

//const char *KEY_DA = "da";
const char *KEY_DA_RAW = "da_raw";

//const char *KEY_DSA_RAW = "dsa_raw";
//const char *KEY_DSA_GC = "dsa_gc";

const char *KEY_BWT_RUNS_FIRST = "bwt_runs_first";
const char *KEY_BWT_RUNS_LAST = "bwt_runs_last";

const char *KEY_TEXT_BWT_RUNS_FIRST = "text_bwt_runs_first";
const char *KEY_TEXT_BWT_RUNS_LAST = "text_bwt_runs_last";

//const char *KEY_R_INDEX = "ri";

//const char *KEY_DOC_END = "doc_end";
const char *KEY_DOC_ISAS = "doc_isas";
//const char *KEY_DOC_DISAS_RAW = "doc_disas_raw";
//const char *KEY_DOC_DISAS_GC = "doc_disas_gc";

//const char *KEY_SADA_RMINQ = "sada_rminq";
//const char *KEY_SADA_RMAXQ = "sada_rmaxq";
//
//const char *KEY_ILCP_BACKWARD = "ilcp_b";
//const char *KEY_ILCP_FORWARD = "ilcp_f";
//const char *KEY_ILCP_BACKWARD_RUN_HEADS = "ilcp_b_run_heads";
//const char *KEY_ILCP_FORWARD_RUN_HEADS = "ilcp_f_run_heads";
//const char *KEY_ILCP_BACKWARD_RMQ = "ilcp_b_rmq";
//const char *KEY_ILCP_FORWARD_RMQ = "ilcp_f_rmq";
//
//const char *KEY_CILCP_BACKWARD_RUN_HEADS = "cilcp_b_run_heads";
//const char *KEY_CILCP_FORWARD_RUN_HEADS = "cilcp_f_run_heads";
//const char *KEY_CILCP_BACKWARD_RMQ = "cilcp_b_rmq";
//const char *KEY_CILCP_FORWARD_RMQ = "cilcp_f_rmq";

const u_int8_t kDocDelimiter = 2;

using BitVector = sdsl::sd_vector<>;
using RangeMinQuery = sdsl::rmq_succinct_sct<true>;
using RangeMaxQuery = sdsl::rmq_succinct_sct<false>;

template<typename II, typename DocBorder, typename DocDelim>
void ConstructDocBorder(II begin, II end, DocBorder &doc_border, const DocDelim &doc_delim, std::size_t size = 0) {
  if (size == 0)
    size = std::distance(begin, end);

  DocBorder tmp_doc_border(size, 0);

  std::size_t i = 0;
  for (auto it = begin; it != end; ++it, ++i) {
    if (*it == doc_delim) {
      tmp_doc_border[i] = 1;
    }
  }

  doc_border.swap(tmp_doc_border);
};

template<uint8_t WIDTH = 8, typename DocBorder, typename DocDelim>
void ConstructDocBorder(const std::string &data_file, DocBorder &doc_border, const DocDelim &doc_delim) {
  sdsl::int_vector_buffer<WIDTH> data_buf(data_file, std::ios::in, 1024 * 1024, WIDTH, true);

  DocBorder tmp_doc_border(data_buf.size(), 0);

  for (std::size_t i = 0; i < data_buf.size(); ++i) {
    if (data_buf[i] == doc_delim) {
      tmp_doc_border[i] = 1;
    }
  }

  doc_border.swap(tmp_doc_border);

  // Previous solution seems to be faster due to iterator comparison.
//  ConstructDocBorder(data_buf.begin(), data_buf.end(), doc_border, doc_delim);
}


//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte- or integer-alphabet and suffix array.
/*!	The algorithm constructs the BWT and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n \log \sigma \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
 *         * conf::KEY_SA
 *  \post BWT exist in the cache. Key
 *         * conf::KEY_BWT for t_width=8 or conf::KEY_BWT_INT for t_width=0
 */
template<uint8_t t_width>
void construct_bwt_and_runs(cache_config &config) {
  static_assert(t_width == 0 or t_width == 8,
                "construct_bwt_and_runs: width must be `0` for integer alphabet and `8` for byte alphabet");

  typedef int_vector<>::size_type size_type;
  typedef int_vector<t_width> text_type;
  typedef int_vector_buffer<t_width> bwt_type;
  const char *KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
  const char *KEY_BWT = key_bwt_trait<t_width>::KEY_BWT;

  //  (1) Load text from disk
  text_type text;
  load_from_cache(text, KEY_TEXT, config);
  size_type n = text.size();
  uint8_t bwt_width = text.width();

  //  (2) Prepare to stream SA from disc and BWT to disc
  size_type buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
  int_vector_buffer<> sa_buf(cache_file_name(conf::KEY_SA, config), std::ios::in, buffer_size);
  std::string bwt_file = cache_file_name(KEY_BWT, config);
  bwt_type bwt_buf(bwt_file, std::ios::out, buffer_size, bwt_width);

  int_vector_buffer<> bwt_runs_first_buf(cache_file_name(KEY_BWT_RUNS_FIRST, config), std::ios::out);
  int_vector_buffer<> text_bwt_runs_first_buf(cache_file_name(KEY_TEXT_BWT_RUNS_FIRST, config), std::ios::out);
  int_vector_buffer<> text_bwt_runs_last_buf(cache_file_name(KEY_TEXT_BWT_RUNS_LAST, config), std::ios::out);

  //  (3) Construct BWT sequentially by streaming SA and random access to text
  size_type to_add[2] = {(size_type) -1, n - 1};
  auto get_text_position = [&sa_buf, &to_add](auto i) {
    auto pos = sa_buf[i];
    return pos + to_add[pos == 0];
  };

  // First BWT value
  bwt_buf[0] = text[get_text_position(0)];
  uint8_t prev_c_bwt = bwt_buf[0];

  size_t nruns = 0; // # BWT runs

  // First position starts the first BWT run.
  bwt_runs_first_buf[nruns] = 0;
  text_bwt_runs_first_buf[nruns] = get_text_position(0);

  auto last_pos = n - 1;
  for (size_type i = 1; i < last_pos; ++i) {
    bwt_buf[i] = text[get_text_position(i)];

    uint8_t c_bwt = bwt_buf[i];
    if (prev_c_bwt != c_bwt) {
      // Last position of the current BWT run
      text_bwt_runs_last_buf[nruns] = get_text_position(i - 1);

      ++nruns;

      // First position of the next BWT run
      bwt_runs_first_buf[nruns] = i;
      text_bwt_runs_first_buf[nruns] = get_text_position(i);

      prev_c_bwt = c_bwt;
    }
  }

  // Last BWT value
  bwt_buf[last_pos] = text[get_text_position(last_pos)];

  // Last position ends the last BWT run
  text_bwt_runs_last_buf[nruns] = get_text_position(last_pos);

  cout << "# BWT runs = " << nruns << endl;

  text_bwt_runs_first_buf.close();
  text_bwt_runs_last_buf.close();
  bwt_runs_first_buf.close();

  bwt_buf.close();

  register_cache_file(KEY_BWT, config);
}

template<uint8_t t_width>
struct sa_trait {
  typedef uint64_t value_type;
  typedef std::vector<value_type> vec_type;
  enum { num_bytes = 0 };
  template<class t_sa>
  static void calc_sa(t_sa &sa, vec_type &text) {
    qsufsort::construct_sa(sa, text);
  }
};

template<>
struct sa_trait<8> {
  typedef uint8_t value_type;
  typedef std::vector<value_type> vec_type;
  enum { num_bytes = 1 };
  template<class t_sa>
  static void calc_sa(t_sa &sa, vec_type &text) {
    algorithm::calculate_sa(text.data(), text.size(), sa);
  }
};

template<uint8_t t_width>
void construct_doc_isa_base(typename sa_trait<t_width>::vec_type &doc_buffer, int_vector<> &doc_isa) {
  int_vector<> sa(doc_buffer.size(), 0, bits::hi(doc_buffer.size()) + 1);
  sa_trait<t_width>::calc_sa(sa, doc_buffer);
  util::bit_compress(sa);
  doc_isa = sa;
  for (std::size_t i = 0; i < doc_buffer.size(); ++i) {
    doc_isa[sa[i]] = i;
  }
}

template<uint8_t t_width>
void construct_doc_isa(cache_config &config, vector<int_vector<>> &doc_isa) {
  typename sa_trait<t_width>::vec_type doc_buffer;
  int_vector_buffer<t_width> text_buf(cache_file_name(conf::KEY_TEXT, config));
  std::size_t doc_id = 0;
  for (std::size_t i = 0; i < text_buf.size(); ++i) {
    if (kDocDelimiter == text_buf[i]) {
      if (doc_buffer.size() > 0) {
        doc_buffer.push_back(0);
        doc_isa.emplace_back(int_vector<>());
        construct_doc_isa_base<t_width>(doc_buffer, doc_isa[doc_isa.size() - 1]);
        ++doc_id;
      }
      doc_buffer.clear();
    } else {
      doc_buffer.push_back(text_buf[i]);
    }
  }
}

template<uint8_t t_width>
auto GetNextDoc(int_vector_buffer<t_width> &_text_buffer,
                typename int_vector_buffer<t_width>::value_type _kDocDelimiter,
                std::size_t &_sp) {
  using ValueType = typename int_vector_buffer<t_width>::value_type;
  std::vector<ValueType> doc_buffer;
  for (; _sp < _text_buffer.size(); ++_sp) {
    ValueType ch = _text_buffer[_sp];

    if (ch == _kDocDelimiter) break;

    doc_buffer.push_back(ch);
  }

  return doc_buffer;
}

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program calculates the SA and BWT for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  construct_config::byte_algo_sa = FLAGS_sais ? SE_SAIS
                                              : LIBDIVSUFSORT; // or LIBDIVSUFSORT for less space-efficient but faster construction

  string data_path = FLAGS_data;

  cache_config config(false, ".", util::basename(FLAGS_data));

  if (!cache_file_exists(conf::KEY_TEXT, config)) {
    int_vector<8> text;
    {
      std::string input;
      {
        std::ifstream fs(data_path);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
      }

      cout << "Text length = " << input.size() << endl;
//        if (any_of(text.begin(), text.end(), [](auto symbol) { return symbol == 0 || symbol == 1; })) {
//            std::cout << "ERROR: Input data contains reserved characters {0x0, 0x1}";
//            exit(1);
//        }

      text.resize(input.size() + 1);

      replace_copy(input.begin(), input.end(), text.begin(), 0, 2);

//        append_zero_symbol(text);
      text[text.size() - 1] = 0;
    }

    store_to_cache(text, conf::KEY_TEXT, config);
    util::clear(text);
  }

  if (!cache_file_exists(conf::KEY_SA, config)) {
    cout << "Calculate Suffix Array ... " << endl;
    construct_sa<8>(config);
    cout << "DONE" << endl;
  }

  if (!cache_file_exists(KEY_DOC_END, config)) {
    std::cout << "Construct Doc Endings ..." << std::endl;

    sdsl::int_vector<8> text;
    load_from_cache(text, conf::KEY_TEXT, config);

    sdsl::bit_vector tmp_doc_endings(text.size(), 0);
    ConstructDocBorder(text.begin(), text.end(), tmp_doc_endings, kDocDelimiter);

    BitVector doc_endings(tmp_doc_endings);
    sdsl::store_to_cache(doc_endings, KEY_DOC_END, config);

    cout << "DONE" << endl;
  }

  if (!cache_file_exists(KEY_DA, config) || !cache_file_exists(KEY_DA_RAW, config)) {
    std::cout << "Construct Document Array (Raw) ..." << std::endl;

    sdsl::int_vector<> sa;
    load_from_cache(sa, conf::KEY_SA, config);

//    sdsl::bit_vector doc_endings;
    BitVector doc_endings;
    load_from_cache(doc_endings, KEY_DOC_END, config);

//    BitVector doc_endings_compact(doc_endings);
//    auto doc_endings_rank = BitVector::rank_1_type(&doc_endings_compact);
    BitVector::rank_1_type doc_endings_rank(&doc_endings);

//    size_t doc_cnt = doc_endings_rank(doc_endings_compact.size());
    size_t doc_cnt = doc_endings_rank(doc_endings.size());

    int_vector<> da(sa.size(), 0, bits::hi(doc_cnt) + 1);
    std::vector<int> da_raw;
    da_raw.reserve(sa.size());
    for (size_t i = 0; i < sa.size(); ++i) {
      da[i] = doc_endings_rank(sa[i]);
      da_raw.emplace_back(da[i]);
    }

    store_to_cache(da, KEY_DA, config);
    {
      auto filepath = cache_file_name(KEY_DA_RAW, config);
      sdsl::osfstream out(filepath, std::ios::binary | std::ios::trunc | std::ios::out);
      serialize_vector(da_raw, out);
    }

    cout << "DONE" << endl;
  }

  if (!cache_file_exists(KEY_DSA_RAW, config)) {
    cout << "Calculate Differential Suffix Array ... " << endl;

    sdsl::int_vector<> sa;
    load_from_cache(sa, conf::KEY_SA, config);

    std::vector<int> sa_diff(sa.size());
    sa_diff[0] = sa[0];
    for (std::size_t i = 1; i < sa.size(); ++i) {
      sa_diff[i] = sa[i] - sa[i - 1];
    }

    auto minmax = std::minmax_element(sa_diff.begin(), sa_diff.end());
    {
      std::ofstream out(cache_file_name(KEY_DSA_RAW, config) + ".info");
      out << sa_diff.size() << std::endl;
      out << *minmax.first << std::endl;
      out << *minmax.second << std::endl;
    }
    std::cout << "Min(SA-Diff): " << *minmax.first << std::endl;
    std::cout << "Max(SA-Diff): " << *minmax.second << std::endl;

    if (*minmax.first < 0) {
      auto minimal = std::abs(*minmax.first);
      std::transform(sa_diff.begin(), sa_diff.end(), sa_diff.begin(), [&minimal](auto v) { return v + minimal; });
    }

    {
      auto mm = std::minmax_element(sa_diff.begin(), sa_diff.end());
      std::cout << "Min(SA-Diff): " << *mm.first << std::endl;
      std::cout << "Max(SA-Diff): " << *mm.second << std::endl;
    }

    {
      auto filepath = cache_file_name(KEY_DSA_RAW, config);
      sdsl::osfstream out(filepath, std::ios::binary | std::ios::trunc | std::ios::out);
      serialize_vector(sa_diff, out);
    }

    cout << "DONE" << endl;
  }

/*
  if (!cache_file_exists(KEY_DSA_GC, config)) {
    cout << "Calculate Differential Suffix Array Grammar Compressed... " << endl;

    grammar::DifferentialSLP<grammar::SLP<>, int_vector<>, int_vector<>, int_vector<>> diff_slp;

    {
      grammar::SLP<> slp;
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(slp);

      std::vector<std::size_t> compact_seq;
      auto report_compact_seq = [&compact_seq](const auto &_var) {
        compact_seq.emplace_back(_var);
      };

      re_pair_reader.Read(cache_file_name(KEY_DSA_RAW, config), slp_wrapper, report_compact_seq);
      std::size_t text_size;
      {
        int_vector_buffer<8> text_buf(cache_file_name(conf::KEY_TEXT, config));
        text_size = text_buf.size();
      }

      uint32_t diff_base;
      {
        std::ifstream in(cache_file_name(KEY_DSA_RAW, config) + ".info");
        int32_t minimal;
        in >> minimal;
        diff_base = minimal < 0 ? std::abs(minimal) : 0;
      }
      std::cout << "Sigma: " << slp.Sigma() << std::endl;
      std::cout << "# Rules: " << slp.GetRules().size() << std::endl;
      std::cout << "Variables: " << slp.Variables() << std::endl;

      auto bit_compress = [](sdsl::int_vector<> &_v) { sdsl::util::bit_compress(_v); };
      diff_slp.Compute(text_size, slp, compact_seq, diff_base, bit_compress, bit_compress, bit_compress);

      auto comp = [&slp](auto it1, auto it2) {
        return (slp.SpanLength(it1) < slp.SpanLength(it2));
      };
      auto mm = std::minmax_element(compact_seq.begin(), compact_seq.end(), comp);
      std::cout << "MinLen: " << *mm.first << " - " << slp.SpanLength(*mm.first) << std::endl;
      std::cout << "MaxLen: " << *mm.second << " - " << slp.SpanLength(*mm.second) << std::endl;
//      std::cout << "Min Root: " << *std::min_element(compact_seq.begin(), compact_seq.end()) << std::endl;
//      for (const auto &var  : compact_seq) {
//        if (var <= diff_slp.Sigma()) {
//          std::cout << "CAGUÉ" << std::endl;
//          break;
//        }
//      }

//      sdsl::int_vector<> sa(text_size, 0, bits::hi(text_size)+1);
//      std::size_t i = 0;
//      auto report = [&sa, &i](auto _suffix) {
//        sa[i] = _suffix;
//        ++i;
//      };
//
//      grammar::ExpandDifferentialSLP(diff_slp, 34962643, 34962643, report);
//
//      sdsl::int_vector<> org_sa;
//      load_from_cache(org_sa, conf::KEY_SA, config);
//
//      std::vector<int> sa_diff(org_sa.size());
//      auto file = cache_file_name(KEY_DSA_RAW, config);
//      sdsl::isfstream in(file, std::ios::binary | std::ios::in);
//      load_vector(sa_diff, in);
//
//      for (std::size_t i = 34962643; i <= 34962643; ++i) {
//        if (org_sa[i] != sa[i]) {
//          std::cout << i << ": " << org_sa[i] << " | " << sa[0] << " | " << sa_diff[i] << " | " << sa_diff[i + 1] << std::endl;
//          break;
//        }
//      }
//      store_to_cache(sa, "sa-test", config);
    }

    sdsl::store_to_cache(diff_slp, KEY_DSA_GC, config);
    cout << "DONE" << endl;
  }
*/

  if (!cache_file_exists(conf::KEY_BWT, config)) {
    cout << "Calculate BWT ... " << endl;
    construct_bwt<8>(config);
//    construct_bwt_and_runs<8>(config);
    cout << "DONE" << endl;
  }

  if (!cache_file_exists(conf::KEY_CSA, config)) {
    cout << "Calculate Compressed Suffix Array ... " << endl;

    sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<63>>, 30, 1000000, sdsl::text_order_sa_sampling<> > csa;
    construct(csa, data_path, config, 8);

    sdsl::store_to_cache(csa, conf::KEY_CSA, config);

    cout << "DONE" << endl;
  }

  if (!cache_file_exists(KEY_DOC_ISAS, config)) {
    cout << "Calculate Inverse Suffix Array for each Doc... " << endl;

    vector<int_vector<> > doc_isa;
    construct_doc_isa<8>(config, doc_isa);

    sdsl::store_to_cache(doc_isa, KEY_DOC_ISAS, config);

    cout << "DONE" << endl;
  }

//  if (!cache_file_exists(KEY_DOC_DISAS_RAW, config) || !cache_file_exists(KEY_DOC_DISAS_GC, config) || true) {
  if (!cache_file_exists(KEY_DOC_DISAS_RAW, config)) {
    cout << "Calculate Raw Inverse Suffix Array for each Doc... " << endl;

    vector<int_vector<> > doc_isas;
    sdsl::load_from_cache(doc_isas, KEY_DOC_ISAS, config);

    std::size_t docs_size = 0;
    for_each(doc_isas.begin(), doc_isas.end(), [&docs_size](const auto &_d) { docs_size += _d.size(); });

    std::vector<int> doc_disas(docs_size);
    int prev = 0;
    std::size_t k = 0;
    for (std::size_t i = 0; i < doc_isas.size(); ++i) {
      for (std::size_t j = 0; j < doc_isas[i].size(); ++j) {
        doc_disas[k++] = doc_isas[i][j] - prev;
        prev = doc_isas[i][j];
      }
    }

    auto minmax = std::minmax_element(doc_disas.begin(), doc_disas.end());
    {
      std::ofstream out(cache_file_name(KEY_DOC_DISAS_RAW, config) + ".info");
      out << doc_disas.size() << std::endl;
      out << *minmax.first << std::endl;
      out << *minmax.second << std::endl;
    }
    std::cout << "Min(SA-Diff): " << *minmax.first << std::endl;
    std::cout << "Max(SA-Diff): " << *minmax.second << std::endl;

    if (*minmax.first < 0) {
      auto min = std::abs(*minmax.first);
      std::transform(doc_disas.begin(), doc_disas.end(), doc_disas.begin(), [&min](auto v) { return v + min; });
    }

    {
      auto mm = std::minmax_element(doc_disas.begin(), doc_disas.end());
      std::cout << "Min(SA-Diff): " << *mm.first << std::endl;
      std::cout << "Max(SA-Diff): " << *mm.second << std::endl;
    }

    {
      auto file = cache_file_name(KEY_DOC_DISAS_RAW, config);
      sdsl::osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
      serialize_vector(doc_disas, out);
    }

    cout << "DONE" << endl;
  }

/*
  if (!cache_file_exists(KEY_DOC_DISAS_GC, config)) {
    cout << "Calculate Differential Inverse Suffix Array Grammar Compressed... " << endl;

    grammar::DifferentialSLP<grammar::SLP<>, int_vector<>, int_vector<>, int_vector<>> diff_slp;
    {
      grammar::SLP<> slp;
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(slp);

      std::vector<std::size_t> compact_seq;
      auto report_compact_seq = [&compact_seq](const auto &_var) {
        compact_seq.emplace_back(_var);
      };

      re_pair_reader.Read(cache_file_name(KEY_DOC_DISAS_RAW, config), slp_wrapper, report_compact_seq);

      std::size_t text_size;
      {
        int_vector_buffer<8> text_buf(cache_file_name(conf::KEY_TEXT, config));
        text_size = text_buf.size();
      }

      uint32_t diff_base;
      {
        std::ifstream in(cache_file_name(KEY_DOC_DISAS_RAW, config) + ".info");
        int32_t minimal;
        in >> minimal;
        diff_base = minimal < 0 ? std::abs(minimal) : 0;
      }
      std::cout << "Sigma: " << slp.Sigma() << std::endl;
      std::cout << "# Rules: " << slp.GetRules().size() << std::endl;
      std::cout << "Variables: " << slp.Variables() << std::endl;

      auto bit_compress = [](sdsl::int_vector<> &_v) { sdsl::util::bit_compress(_v); };
      diff_slp.Compute(text_size - 1, slp, compact_seq, diff_base, bit_compress, bit_compress, bit_compress);

      vector<int_vector<> > doc_isas;
      load_from_cache(doc_isas, KEY_DOC_ISAS, config);

      vector<int_vector<> > doc_isas_1;
      doc_isas_1.reserve(doc_isas.size());

      std::size_t k = 0;
      std::size_t i = 0;
      auto report = [&doc_isas_1, &i, &k](auto _v) {
        doc_isas_1[doc_isas_1.size() - 1][i] = _v;
        ++i;
        ++k;
      };

      for (std::size_t j = 0; j < doc_isas.size(); ++j) {
        i = 0;
        doc_isas_1.emplace_back(int_vector<>(doc_isas[j].size()));
        grammar::ExpandDifferentialSLP(diff_slp, k, k + doc_isas[j].size() - 1, report);
        util::bit_compress(doc_isas_1[doc_isas_1.size() - 1]);
      }

      store_to_cache(doc_isas_1, "test", config);
    }

    sdsl::store_to_cache(diff_slp, KEY_DOC_DISAS_GC, config);
    cout << "DONE" << std::endl;
  }

//  }

  if (!cache_file_exists(KEY_ILCP_BACKWARD, config) || !cache_file_exists(KEY_ILCP_FORWARD, config)
      || !cache_file_exists(KEY_ILCP_BACKWARD_RMQ, config) || !cache_file_exists(KEY_ILCP_FORWARD_RMQ, config)
      || !cache_file_exists(KEY_CILCP_BACKWARD_RMQ, config) || !cache_file_exists(KEY_CILCP_FORWARD_RMQ, config)) {
    int_vector<> ilcps[2]; // [0] == backward && [1] == forward

    if (!load_from_cache(ilcps[0], KEY_ILCP_BACKWARD, config) || !load_from_cache(ilcps[1], KEY_ILCP_FORWARD, config)) {
      cout << "Calculate ILCP... " << endl;

      std::size_t text_size;
      std::vector<int_vector<>> lcp_docs;
      {
        int_vector_buffer<8> text_buf(cache_file_name(conf::KEY_TEXT, config));
        text_size = text_buf.size();

//        // Collection ISA
//        sdsl::construct_isa(config);
//        int_vector_buffer<> isa_buf(cache_file_name(conf::KEY_ISA, config));

        cache_config config_doc(false, ".", util::basename(FLAGS_data) + "-doc");

        std::size_t pos = 0;
        std::size_t doc = 0;
        while (pos < text_buf.size()) {
//          std::cout << "Doc " << doc++ << std::endl;

          // Document text
//          std::cout << "  Text" << std::endl;
//          std::size_t sp_doc = pos;
          {
            sdsl::int_vector<8> text_doc;
            {
              auto doc_buffer = GetNextDoc(text_buf, 2, pos);
              ++pos; // Set at first position of next document

              text_doc.resize(doc_buffer.size() + 1);

              std::copy(doc_buffer.begin(), doc_buffer.end(), text_doc.begin());
              text_doc[text_doc.size() - 1] = 0;
            }
            sdsl::store_to_cache(text_doc, conf::KEY_TEXT, config_doc);
          }

          // Document SA
//          std::cout << "  SA" << std::endl;
          sdsl::construct_sa<8>(config_doc);

          // Document ISA
//          std::cout << "  ISA" << std::endl;
          sdsl::construct_isa(config_doc);

          // Document LCP
//          std::cout << "  LCP" << std::endl;
          sdsl::construct_lcp_kasai<8>(config_doc);
          lcp_docs.emplace_back(sdsl::int_vector<>());
          sdsl::load_from_cache(lcp_docs[lcp_docs.size() - 1], conf::KEY_LCP, config_doc);

//        // Construct ILCP
//        std::cout << "  ILCP" << std::endl;
//        sdsl::int_vector_buffer<> lcp_doc_buf(cache_file_name(conf::KEY_LCP, config_doc));
//        sdsl::int_vector_buffer<> sa_doc_buf(cache_file_name(conf::KEY_SA, config_doc));
//        isa_buf[sp_doc]; // Set isa_buf at first suffix of current document.
//        for (std::size_t i = 0; i < sa_doc_buf.size(); ++i) {
//          ilcp_backward[isa_buf[sa_doc_buf[i] + sp_doc]] = lcp_doc_buf[i];
//        }

          for (const auto &item : config_doc.file_map) {
//            std::cout << item.first << "---" << item.second << std::endl;
            std::remove(item.second.c_str());
          }
        }
      }

      // Construct ILCP
      std::cout << "  ILCP" << std::endl;

      ilcps[0].resize(text_size);
      ilcps[1].resize(text_size);
      sdsl::int_vector_buffer<> da(cache_file_name(KEY_DA, config));
      std::vector<std::size_t> doc_suffix_indices(lcp_docs.size(), 0);

      for (std::size_t i = 0; i < da.size(); ++i) {
        std::size_t doc = da[i];
        auto &doc_suffix_idx = doc_suffix_indices[doc];
        const auto &doc_lcp = lcp_docs[doc];

        ilcps[0][i] = doc_lcp[doc_suffix_idx];
        ilcps[1][i] = doc_lcp[++doc_suffix_idx % doc_lcp.size()];
      }

      sdsl::store_to_cache(ilcps[0], KEY_ILCP_BACKWARD, config);
      sdsl::store_to_cache(ilcps[1], KEY_ILCP_FORWARD, config);

      cout << "DONE" << endl;
    }

    if (!cache_file_exists(KEY_ILCP_BACKWARD_RMQ, config) || !cache_file_exists(KEY_ILCP_FORWARD_RMQ, config)) {
      cout << "Calculate ILCP Left/Right RMQ... " << endl;

      const char *key_ilcps_run_beginnings[2] = {KEY_ILCP_BACKWARD_RUN_HEADS, KEY_ILCP_FORWARD_RUN_HEADS};
      const char *key_ilcps_rmq[2] = {KEY_ILCP_BACKWARD_RMQ, KEY_ILCP_FORWARD_RMQ};

      for (int direction = 0; direction < 2; ++direction) {
        std::cout << "Direction: " << direction << std::endl;

        const auto &ilcp = ilcps[direction];

        sdsl::bit_vector ilcp_run_beginnings(ilcp.size(), 0);
        std::vector<std::size_t> ilcp_runs;

        ilcp_runs.emplace_back(ilcp[0]);
        ilcp_run_beginnings[0] = 1;
        for (std::size_t i = 1; i < ilcp.size(); ++i) {
          if (ilcp[i - 1] != ilcp[i]) {
            ilcp_runs.emplace_back(ilcp[i]);
            ilcp_run_beginnings[i] = 1;
          }
        }

//        {
//          std::ofstream out(std::to_string(direction) + "-tmp.txt");
//          for (int i =31892292; i <= 31892408; ++i) {
//            out << i << " - " << ilcp_runs[i] << std::endl;
//          }
//        }

        store_to_cache(ilcp_run_beginnings, key_ilcps_run_beginnings[direction], config);

        std::cout << "# ILCP runs: " << ilcp_runs.size() << std::endl;

        RangeMinQuery range_min_query(&ilcp_runs);
        store_to_cache(range_min_query, key_ilcps_rmq[direction], config);
      }

      cout << "DONE" << endl;
    }

    if (!cache_file_exists(KEY_CILCP_BACKWARD_RMQ, config) || !cache_file_exists(KEY_CILCP_FORWARD_RMQ, config)) {
      cout << "Calculate CILCP Left/Right RMQ... " << endl;

      const char *key_ilcps_run_heads[2] = {KEY_CILCP_BACKWARD_RUN_HEADS, KEY_CILCP_FORWARD_RUN_HEADS};
      const char *key_ilcps_rmq[2] = {KEY_CILCP_BACKWARD_RMQ, KEY_CILCP_FORWARD_RMQ};

      sdsl::int_vector<> da;
      sdsl::load_from_cache(da, KEY_DA, config);

//      {
//        std::ofstream out("tmp-ilcp.txt");
//        for (int i = 36854822; i <= 36856834; ++i) {
//          out << i << " - " << da[i] << " - " << ilcps[1][i] << std::endl;
//        }
//      }

      for (int direction = 0; direction < 2; ++direction) {
        std::cout << "Direction: " << direction << std::endl;
        const auto &ilcp = ilcps[direction];

        sdsl::bit_vector ilcp_run_heads(ilcp.size(), 0);
        std::vector<std::size_t> ilcp_runs;

        ilcp_runs.emplace_back(ilcp[0]);
        ilcp_run_heads[0] = 1;
        for (std::size_t i = 1; i < ilcp.size(); ++i) {
          std::size_t l = ilcp[i - 1];
          std::size_t d = da[i - 1];
          if (l == ilcp[i]) {
            while (++i < ilcp.size() && l == ilcp[i]) {}
          } else if (d == da[i]) {
            do {
              l = std::min(l, ilcp[i]);
            } while (++i < ilcp.size() && d == da[i]);
          }

          if (i < ilcp.size()) {
            ilcp_runs[ilcp_runs.size() - 1] = l;
            ilcp_runs.emplace_back(ilcp[i]);
            ilcp_run_heads[i] = 1;
          }
        }

//        {
//          sdsl::bit_vector::select_1_type select(&ilcp_run_heads);
//
//          std::ofstream out(std::to_string(direction) + "-tmp-cilcp-runs.txt");
//          for (int i = 904788; i <= 904796; ++i) {
//            out << i << " - " << select(i + 1) << " - " << ilcp_runs[i] << std::endl;
//          }
//        }

        store_to_cache(ilcp_run_heads, key_ilcps_run_heads[direction], config);

        std::cout << "# CILCP runs: " << ilcp_runs.size() << std::endl;

        RangeMinQuery range_min_query(&ilcp_runs);
        store_to_cache(range_min_query, key_ilcps_rmq[direction], config);
      }

      cout << "DONE" << endl;
    }
  }
*/

  if (!cache_file_exists(KEY_R_INDEX, config)) {
    std::cout << "Construct RI ..." << std::endl;

    std::string input;
    {
      std::ifstream fs(data_path);
      std::stringstream buffer;
      buffer << fs.rdbuf();

      input = buffer.str();
    }

    std::replace(input.begin(), input.end(), '\0', '\2');
    ri::r_index<> r_idx(input, FLAGS_sais);

    sdsl::store_to_cache(r_idx, KEY_R_INDEX, config);

    cout << "DONE" << endl;
  }

/*
  if (!cache_file_exists(KEY_SADA_RMINQ, config) || !cache_file_exists(KEY_SADA_RMAXQ, config)) {
    std::cout << "Construct SADA Range Min/Max Query ..." << std::endl;

    sdsl::int_vector<> da;
    load_from_cache(da, KEY_DA, config);

    size_t doc_cnt;
    {
      sdsl::bit_vector doc_endings;
      load_from_cache(doc_endings, KEY_DOC_END, config);
      using BitVector = sdsl::sd_vector<>;

      BitVector doc_endings_compact(doc_endings);
      auto doc_endings_rank = BitVector::rank_1_type(&doc_endings_compact);

      doc_cnt = doc_endings_rank(doc_endings_compact.size());
    }

    // Compute range minimum query on previous document array
    {
      sdsl::int_vector<> prev_docs;
      dret::ConstructPrevDocArray(da, doc_cnt, prev_docs);
      RangeMinQuery range_min_query(&prev_docs);
      store_to_cache(range_min_query, KEY_SADA_RMINQ, config);
    }

    // Compute range maximum query on next document array
    {
      sdsl::int_vector<> next_docs;
      dret::ConstructNextDocArray(da, doc_cnt, next_docs);
      RangeMaxQuery range_max_query(&next_docs);
      store_to_cache(range_max_query, KEY_SADA_RMAXQ, config);
    }
    cout << "DONE" << endl;
  }
*/
  return 0;
}
