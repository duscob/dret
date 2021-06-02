//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/5/20.
//

#include <iostream>
#include <fstream>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>

#include "../benchmark/factory.h"
#include "dret/doc_freq_index_sada.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_name, "data", "Data file basename.");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

using BitVectorCompact = sdsl::sd_vector<>;
using BitVectorCompactRank = BitVectorCompact::rank_1_type;
using BitVectorCompactSelect = BitVectorCompact::select_1_type;

auto compute_sa_chunk_by_doc(const sdsl::int_vector<> &t_da, const std::pair<ulint, ulint> t_range) {
  std::map<uint16_t, std::vector<std::size_t>> pos_by_doc;
  for (auto i = t_range.first; i <= t_range.second; ++i) {
    auto d = t_da[i];
    if (pos_by_doc.find(d) != pos_by_doc.end())
      continue;

    {
      std::size_t c = 3;
      auto j = t_range.first - 1;
      while (0 < c) {
        while (t_da[j] != d) { --j; }

        pos_by_doc[d].emplace_back(j);

        --j;
        --c;
      }
    }

    for (auto j = i; j <= t_range.second; ++j) {
      if (t_da[j] != d)
        continue;

      pos_by_doc[d].emplace_back(j);
    }

    {
      std::size_t c = 3;
      auto j = t_range.second + 1;
      while (0 < c) {
        while (t_da[j] != d) { ++j; }

        pos_by_doc[d].emplace_back(j);

        ++j;
        --c;
      }
    }
  }

  return pos_by_doc;
}

auto compute_da_grammar_tree(const grammar::SLP<> &t_slp,
                             const std::vector<std::size_t> &t_slp_roots,
                             std::size_t t_sp,
                             std::size_t t_ep) {
  std::cout << "No Rules: " << t_slp.Start() << std::endl;
  std::cout << "No Roots: " << t_slp_roots.size() << std::endl;
  std::cout << "[" << t_sp << ":" << t_ep << "]" << std::endl;

  std::size_t root_offset = t_sp;
  std::size_t root_idx = 0;
  std::size_t prev_space = 0;

  while (root_idx < t_slp_roots.size() && t_slp.SpanLength(t_slp_roots[root_idx]) <= root_offset) {
    prev_space += t_slp.SpanLength(t_slp_roots[root_idx]);
    root_offset -= t_slp.SpanLength(t_slp_roots[root_idx]);
    ++root_idx;
  }

  auto root = t_slp_roots[root_idx];
  std::cout << root_idx << " : " << root << "[" << t_slp.SpanLength(root) << "] - " << root_offset << " - "
            << prev_space << std::endl;
}

auto print_padding(std::ostream &t_out, std::size_t t_n) {
  for (int i = 0; i < t_n; ++i) {
    t_out << " ";
  }
}

void print_grammar_tree(const grammar::SLP<> &t_slp,
                        std::size_t t_v,
                        std::size_t t_level,
                        const std::string &t_id,
                        std::size_t &t_n_leaves,
                        std::ostream &t_out,
                        std::ostream &t_out_styles) {
  print_padding(t_out, t_level);

  auto id = t_id + std::to_string(t_v);

  if (t_slp.IsTerminal(t_v)) {
//    t_out << t_v << "/" << id << std::endl;
//    t_out << t_v << std::endl;
    t_out << "XXX" << t_n_leaves << std::endl;
    ++t_n_leaves;
    return;
  }

  t_out << "\"" << id << "\"/" << t_v << " [" << id << "sty]" << " -- {" << std::endl;
//  t_out << t_v << "[fresh nodes] -- {" << std::endl;

  t_out_styles << "\\tikzset{" << id << "sty/.style={}}" << std::endl;

  auto children = t_slp[t_v];
  print_grammar_tree(t_slp, children.first, t_level + 1, id + "l", t_n_leaves, t_out, t_out_styles);

  print_padding(t_out, t_level + 1);
  t_out << "," << std::endl;

  print_grammar_tree(t_slp, children.second, t_level + 1, id + "r", t_n_leaves, t_out, t_out_styles);

  print_padding(t_out, t_level);
  t_out << "}" << std::endl;
}

auto compute_da_grammar_tree(const sdsl::int_vector<> &t_da,
                             std::pair<std::size_t, std::size_t> t_range,
                             std::ostream &t_out,
                             std::ostream &t_out_styles,
                             std::ostream &t_csv) {
  auto bit = t_da.begin() + t_range.first;
  auto eit = t_da.begin() + t_range.second + 1;

  grammar::SLP<> slp;
  {
    grammar::RePairEncoder<true> grammar_encoder;
    auto slp_wrapper = grammar::BuildSLPWrapper(slp);

    grammar_encoder.Encode(bit, eit, slp_wrapper);
  }

  auto seq = slp.Span(slp.Start());

  auto print_seq = [](auto bit, auto eit) {
    for (auto it = bit; it != eit; ++it) {
      std::cout << *it << " ";
    }
    std::cout << std::endl;
  };

  print_seq(bit, eit);
  print_seq(seq.begin(), seq.end());

  std::size_t n_leaves = 0;
  print_grammar_tree(slp, slp.Start(), 0, "", n_leaves, t_out, t_out_styles);

  t_csv << "Var,BV_F,BV_L" << std::endl;
//  sdsl::bit_vector bv_first;
//  sdsl::bit_vector bv_last;
  auto values = dret::BuildFirstAndLastOccs<sdsl::bit_vector>(slp);
  t_csv << 14 << ",\"" << values.first[14] << "\",\"" << values.second[14] << "\"" << std::endl;

}

int main(int argc, char *argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_patterns.empty() || FLAGS_data_name.empty() || FLAGS_data_dir.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::cache_config config(true, FLAGS_data_dir, FLAGS_data_name);

  ri::r_index<> r_idx;
  Load(r_idx, KEY_R_INDEX, config, "R-Index");

  sdsl::int_vector<8> text;
  Load(text, sdsl::conf::KEY_TEXT, config, "Text");

  sdsl::int_vector<> sa;
  Load(sa, sdsl::conf::KEY_SA, config, "SA");

  sdsl::int_vector<> da;
  Load(da, KEY_DA, config, "DA");
//
////  std::vector<int> seq = {0, 1, 0, 1, 0, 1};
//  std::vector<int> seq = {2, 1, 2, 1, 2, 1};
//
//  auto bit = da.begin();
//  auto it = bit;
//
//  while ((it = std::search(bit, da.end(), seq.begin(), seq.end())) != da.end()) {
//    auto os = it - da.begin();
//    if (it != da.end() && std::count(text.begin() + sa[os], text.begin() + sa[os] + 20, '-') == 0) {
//      std::cout << "Y --- " << it - da.begin() << std::endl;
//      std::vector<int> seq = {0, 1, 0, 1, 0, 1};
//      for (int i = -10; i < 10; ++i) {
//        auto j = os + i;
//        std::cout << j << "\t" << da[j] << "\t" << sa[j] << "\t|";
//        std::copy_n(text.begin() + sa[j], 20, std::ostream_iterator<uint8_t>(std::cout, ""));
//        std::cout << std::endl;
//      }
//
////      break;
//    } else {
//      std::cout << "N" << std::endl;
//    }
//
//    bit = it + 1;
//  }

  BitVectorCompact doc_endings;
  Load(doc_endings, KEY_DOC_END, config, "Document Endings");
  auto doc_endings_rank = BitVectorCompactRank(&doc_endings);
  auto doc_endings_select = BitVectorCompactSelect(&doc_endings);

//  auto doc_cnt = doc_endings_rank(doc_endings.size());
//  // Compute range minimum query on previous document array
//  {
//    sdsl::int_vector<> prev_docs;
//    dret::ConstructPrevDocArray(da, doc_cnt, prev_docs);
//    store_to_cache(prev_docs, "sada_c_prev", config);
//  }

//  // Compute range maximum query on next document array
//  {
//    sdsl::int_vector<> next_docs;
//    dret::ConstructNextDocArray(da, doc_cnt, next_docs);
//    store_to_cache(next_docs, "sada_c_next", config);
//  }

  sdsl::int_vector<> sada_c[2]; // [0] == backward && [1] == forward
  Load(sada_c[0], "sada_c_prev", config, "SADA C Backward");
//  Load(sada_c[1], "sada_c_next", config, "SADA C Forward");

  sdsl::int_vector<> ilcps[2]; // [0] == backward && [1] == forward
  Load(ilcps[0], KEY_ILCP_BACKWARD, config, "ILCP Backward");
//  Load(ilcps[1], KEY_ILCP_FORWARD, config, "ILCP Forward");
  sdsl::bit_vector ilcp_run_heads;
  Load(ilcp_run_heads, KEY_ILCP_BACKWARD_RUN_HEADS, config, "ILCP Run Heads Backward");
  sdsl::bit_vector cilcp_run_heads;
  Load(cilcp_run_heads, KEY_CILCP_BACKWARD_RUN_HEADS, config, "CILCP Run Heads Backward");

////  auto it = std::find(ilcp_run_heads.begin() + 500000, ilcp_run_heads.end(), 0);
//  std::vector<int> seq = {0, 0};
////  std::vector<int> seq = {0, 1, 1, 1, 0};
////  std::vector<int> seq = {0, 1, 1, 0};
////  std::vector<int> seq = {0, 1, 0, 1, 1, 1, 1, 1, 0};
////  std::vector<int> seq = {0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0};
//  auto it = std::search(ilcp_run_heads.begin() + 800000, ilcp_run_heads.end(), seq.begin(), seq.end());
//  if (it != ilcp_run_heads.end()) {
//    auto os = it - ilcp_run_heads.begin();
//    std::cout << "Y --- " << os << std::endl;
//    for (int i = -10; i < 10; ++i) {
//      auto j = os + i;
//      std::cout << j << "\t" << da[j] << "\t" << sa[j] << "\t" << ilcps[0][j] << "\t" << ilcp_run_heads[j] << "\t";
//      std::copy_n(text.begin() + sa[j], 20, std::ostream_iterator<uint8_t>(std::cout, ""));
//      std::cout << "\t";
//      std::copy_n(text.begin() + sa[j], 15, std::ostream_iterator<int>(std::cout, " "));
//      std::cout << std::endl;
//    }
//  } else {
//    std::cout << "N --- " << std::endl;
//  }

  sdsl::wt_ap<> wt_da_ap;
  Load(wt_da_ap, KEY_DA_WT + std::string("_ap"), config, "WT DA");

  grammar::SLP<> da_slp;
  std::vector<std::size_t> da_slp_roots;
  {
    grammar::RePairReader<false> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(da_slp);

    auto report_compact_seq = [&da_slp_roots](const auto &_var) {
      da_slp_roots.emplace_back(_var);
    };

    re_pair_reader.Read(cache_file_name(KEY_DA_RAW, config), slp_wrapper, report_compact_seq);
  }

  // Query patterns
  std::vector<std::string> patterns;
  {
    std::ifstream pattern_file(FLAGS_patterns.c_str(), std::ios_base::binary);
    if (!pattern_file) {
      std::cerr << "ERROR: Failed to open patterns file!" << std::endl;
      return 3;
    }

    std::string pat;
    while (std::getline(pattern_file, pat)) {
      if (pat.empty())
        continue;

      patterns.emplace_back(pat);
      std::ofstream file_info(pat + ".info.csv");
      file_info << "pattern,range_begin,range_end,erange_begin,erange_end" << std::endl;
      file_info << pat;
      auto range = r_idx.count(pat);
//      file_info << range.first << " " << range.second << std::endl;
//      file_info << range.second - range.first + 1 << std::endl;
      file_info << "," << range.first << "," << range.second;

      auto c = range.second - range.first + 1;
      auto size = 26ul;
      auto offset = c < size ? (size - c) / 2 + 1 : 5;

      auto range_e = std::make_pair(std::max(range.first, offset) - offset,
                                    std::min(range.second + offset, r_idx.text_size() - 1));
//      file_info << range_e.first << " " << range_e.second << std::endl;
//      file_info << range_e.second - range_e.first + 1 << std::endl;
      file_info << "," << range_e.first << "," << range_e.second;

      {
        std::ofstream file_data(pat + ".csv");
        file_data << "Position,SA,Suffix,DA,SADA C Backward,ILCP Backward,ILCP Run Heads,CILCP Run Heads" << std::endl;
        for (auto i = range_e.first; i <= range_e.second; ++i) {
          file_data << i;
          auto sa_value = sa[i];
          file_data << "," << sa[i];
          file_data << ",";
          std::copy_n(text.begin() + sa_value, 15, std::ostream_iterator<uint8_t>(file_data, ""));
          file_data << "";
//        file_data << "," << da[i];
          file_data << "," << doc_endings_rank(sa_value);
          file_data << "," << sada_c[0][i];
//        file_data << "," << sada_c[1][i];
          file_data << "," << ilcps[0][i];
//        file_data << "," << ilcps[1][i];
          file_data << "," << ilcp_run_heads[i];
          file_data << "," << cilcp_run_heads[i];

          file_data << std::endl;
        }
      }

      {
        auto pos_by_doc = compute_sa_chunk_by_doc(da, range);
        for (const auto &item: pos_by_doc) {
          uint16_t doc = item.first;
//          std::cout << doc << std::endl;
          std::ofstream fileout(pat + "." + std::to_string(doc) + ".csv");

          fileout << "Position,Position_d,SA,SA_d" << std::endl;

          auto pos_d = wt_da_ap.rank(item.second.front(), item.first);
          auto start_pos_doc = 0 < doc ? doc_endings_select(doc) : 0;
//          std::cout << start_pos_doc << std::endl;
          for (const auto &pos : item.second) {
            auto sa_value = sa[pos];
            fileout << pos << "," << ++pos_d << "," << sa_value << "," << sa_value - start_pos_doc;
            fileout << std::endl;
//            std::cout << pos << "," << sa_value << "," << sa_value - start_pos_doc;
          }
        }
      }


      // Grammar Tree DA
//      compute_da_grammar_tree(da_slp, da_slp_roots, range.first, range.second);
      {
        std::ofstream fileout(pat + ".grm");
        std::ofstream fileout_sty(pat + ".grm.sty");
        std::ofstream fileout_csv(pat + ".grm.csv");
        compute_da_grammar_tree(da, range_e, fileout, fileout_sty, fileout_csv);
      }
    }

    pattern_file.close();
  }

  return 0;
}