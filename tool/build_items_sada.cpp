//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/11/20.
//

#include <iostream>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>
#include <sdsl/util.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/rmq_support.hpp>

#include <dret/doc_freq_index_sada.h>

#include "definitions.h"


DEFINE_string(data, "", "Data file. (MANDATORY)");

using RangeMinQuery = sdsl::rmq_succinct_sct<true>;
using RangeMaxQuery = sdsl::rmq_succinct_sct<false>;

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program calculates the SA and BWT for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

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
    std::cout << "DONE" << std::endl;
  }

  return 0;
}
