//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/11/20.
//

#include <iostream>

#include <gflags/gflags.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/construct_isa.hpp>
#include <sdsl/construct_lcp.hpp>

#include <sdsl/rmq_support.hpp>

#include "definitions.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");

using BitVector = sdsl::sd_vector<>;
using RangeMinQuery = sdsl::rmq_succinct_sct<true>;

template<uint8_t t_width>
auto GetNextDoc(sdsl::int_vector_buffer<t_width> &_text_buffer,
                typename sdsl::int_vector_buffer<t_width>::value_type _kDocDelimiter,
                std::size_t &_sp) {
  using ValueType = typename sdsl::int_vector_buffer<t_width>::value_type;
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

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

  sdsl::int_vector<> ilcps[2]; // [0] == backward && [1] == forward

  if (!load_from_cache(ilcps[0], KEY_ILCP_BACKWARD, config) || !load_from_cache(ilcps[1], KEY_ILCP_FORWARD, config)) {
    std::cout << "Calculate ILCP... " << std::endl;

    std::size_t text_size;
    std::vector<sdsl::int_vector<>> lcp_docs;
    {
      std::cout << "  Document LCPs" << std::endl;
      sdsl::int_vector_buffer<8> text_buf(cache_file_name(sdsl::conf::KEY_TEXT, config));
      text_size = text_buf.size();

      sdsl::cache_config config_doc(false, ".", sdsl::util::basename(FLAGS_data) + "-doc");

      std::size_t pos = 0;
      while (pos < text_buf.size()) {
        // Document text
        {
          sdsl::int_vector<8> text_doc;
          {
            auto doc_buffer = GetNextDoc(text_buf, 2, pos);
            ++pos; // Set at first position of next document

            text_doc.resize(doc_buffer.size() + 1);

            std::copy(doc_buffer.begin(), doc_buffer.end(), text_doc.begin());
            text_doc[text_doc.size() - 1] = 0;
          }
          sdsl::store_to_cache(text_doc, sdsl::conf::KEY_TEXT, config_doc);
        }

        // Document SA
        sdsl::construct_sa<8>(config_doc);

        // Document ISA
        sdsl::construct_isa(config_doc);

        // Document LCP
        sdsl::construct_lcp_kasai<8>(config_doc);
        lcp_docs.emplace_back(sdsl::int_vector<>());
        sdsl::load_from_cache(lcp_docs[lcp_docs.size() - 1], sdsl::conf::KEY_LCP, config_doc);

        for (const auto &item : config_doc.file_map) {
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

    std::cout << "DONE" << std::endl;
  }


  if (!cache_file_exists(KEY_ILCP_BACKWARD_RMQ, config) || !cache_file_exists(KEY_ILCP_FORWARD_RMQ, config)) {
    std::cout << "Calculate ILCP Left/Right RMQ... " << std::endl;

    const char *key_ilcps_run_heads[2] = {KEY_ILCP_BACKWARD_RUN_HEADS, KEY_ILCP_FORWARD_RUN_HEADS};
    const char *key_ilcps_rmq[2] = {KEY_ILCP_BACKWARD_RMQ, KEY_ILCP_FORWARD_RMQ};

    for (int direction = 0; direction < 2; ++direction) {
      std::cout << "  Direction: " << direction << std::endl;

      const auto &ilcp = ilcps[direction];

      sdsl::bit_vector ilcp_run_heads(ilcp.size(), 0);
      std::vector<std::size_t> ilcp_runs;

      ilcp_runs.emplace_back(ilcp[0]);
      ilcp_run_heads[0] = 1;
      for (std::size_t i = 1; i < ilcp.size(); ++i) {
        if (ilcp[i - 1] != ilcp[i]) {
          ilcp_runs.emplace_back(ilcp[i]);
          ilcp_run_heads[i] = 1;
        }
      }

      std::cout << "    # ILCP runs: " << ilcp_runs.size() << std::endl;
      store_to_cache(ilcp_run_heads, key_ilcps_run_heads[direction], config);

      RangeMinQuery range_min_query(&ilcp_runs);
      store_to_cache(range_min_query, key_ilcps_rmq[direction], config);
    }

    std::cout << "DONE" << std::endl;
  }


  if (!cache_file_exists(KEY_CILCP_BACKWARD_RMQ, config) || !cache_file_exists(KEY_CILCP_FORWARD_RMQ, config)) {
    std::cout << "Calculate CILCP Left/Right RMQ... " << std::endl;

    const char *key_ilcps_run_heads[2] = {KEY_CILCP_BACKWARD_RUN_HEADS, KEY_CILCP_FORWARD_RUN_HEADS};
    const char *key_ilcps_rmq[2] = {KEY_CILCP_BACKWARD_RMQ, KEY_CILCP_FORWARD_RMQ};

    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, KEY_DA, config);

    for (int direction = 0; direction < 2; ++direction) {
      std::cout << "  Direction: " << direction << std::endl;
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

      std::cout << "    # CILCP runs: " << ilcp_runs.size() << std::endl;
      store_to_cache(ilcp_run_heads, key_ilcps_run_heads[direction], config);

      RangeMinQuery range_min_query(&ilcp_runs);
      store_to_cache(range_min_query, key_ilcps_rmq[direction], config);
    }

    std::cout << "DONE" << std::endl;
  }

  return 0;
}
