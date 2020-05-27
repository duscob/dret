//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/26/20.
//

#include <iostream>

#include <boost/filesystem.hpp>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>
#include <sdsl/bit_vectors.hpp>

#include <grammar/slp.h>
#include <grammar/re_pair.h>
#include <grammar/slp_helper.h>

#include <dret/tf.h>

#include "definitions.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program build items of GCDA index.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));
  auto datafile = cache_file_name(KEY_DA_RAW, config);

  if (!sdsl::cache_file_exists(KEY_GCDA_FIRST_OCCS, config)
      || !sdsl::cache_file_exists(KEY_GCDA_LAST_OCCS, config)
      || !sdsl::cache_file_exists(KEY_GCDA_CSEQ_HEADS, config)) {

    grammar::SLP<> slp;
    std::vector<std::size_t> compact_seq;
    {
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(slp);

      auto report_compact_seq = [&compact_seq](const auto &_var) {
        compact_seq.emplace_back(_var);
      };

      re_pair_reader.Read(datafile, slp_wrapper, report_compact_seq);
    }

    if (!sdsl::cache_file_exists(KEY_GCDA_FIRST_OCCS, config) || !sdsl::cache_file_exists(KEY_GCDA_LAST_OCCS, config)) {
      using OccsBV = std::vector<sdsl::bit_vector>;
      auto occs = std::make_pair(OccsBV(slp.Variables() + 1), OccsBV(slp.Variables() + 1));

      for (auto i = slp.Start(); slp.Sigma() < i; --i) {
        dret::BuildFirstAndLastOccs(slp, i, occs.first, occs.second);
      }

      sdsl::store_to_cache(occs.first, KEY_GCDA_FIRST_OCCS, config);
      sdsl::store_to_cache(occs.second, KEY_GCDA_LAST_OCCS, config);
    }

    if (!sdsl::cache_file_exists(KEY_GCDA_CSEQ_HEADS, config)) {
      std::size_t seq_size = boost::filesystem::file_size(datafile) / 4;

      sdsl::bit_vector cseq_heads(seq_size, 0);
      std::size_t pos = 0;
      for (const auto &var : compact_seq) {
        cseq_heads[pos] = 1;
        pos += slp.SpanLength(var);
      }

      sdsl::store_to_cache(cseq_heads, KEY_GCDA_CSEQ_HEADS, config);
    }
  }

  return 0;
}