//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/26/20.
//

#include <iostream>
#include <sys/stat.h>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>
#include <sdsl/bit_vectors.hpp>

#include <grammar/slp.h>
#include <grammar/re_pair.h>
#include <grammar/slp_helper.h>
#include <grammar/sampled_slp.h>

#include <dret/tf.h>

#include "definitions.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_int32(bs, 512, "Block size.");
DEFINE_int32(sf, 4, "Storing factor.");


auto GetFileSize(const std::string &filename)
{
  struct stat64 stat_buf;
  int rc = stat64(filename.c_str(), &stat_buf);
  return rc == 0 ? stat_buf.st_size : -1;
}


int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program build items of GCDA index.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

  auto bit_compress = [](sdsl::int_vector<> &_v) { sdsl::util::bit_compress(_v); };

  if (!sdsl::cache_file_exists(KEY_GCDA_FIRST_OCCS, config)
      || !sdsl::cache_file_exists(KEY_GCDA_LAST_OCCS, config)
      || !sdsl::cache_file_exists(KEY_GCDA_CSEQ_HEADS, config)) {

    auto datafile = cache_file_name(KEY_DA_RAW, config);
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
      std::size_t seq_size = GetFileSize(datafile) / 4;

      sdsl::bit_vector cseq_heads(seq_size, 0);
      std::size_t pos = 0;
      for (const auto &var : compact_seq) {
        cseq_heads[pos] = 1;
        pos += slp.SpanLength(var);
      }

      sdsl::store_to_cache(cseq_heads, KEY_GCDA_CSEQ_HEADS, config);
    }
  }

  grammar::CombinedSLP<> cslp;
  if (!sdsl::load_from_cache(cslp, KEY_GCDA_CSLP, config)
      || !sdsl::cache_file_exists(KEY_GCDA_CSLP_DOCS, config)
      || !sdsl::cache_file_exists(KEY_GCDA_CSLP_OCCS, config)) {
    std::cout << "Construct Combined SLP of Document Array (Raw) ..." << std::endl;
    auto datafile = cache_file_name(KEY_DA_RAW, config);
    grammar::SLP<> slp;
    std::vector<std::size_t> compact_seq;
    {
      grammar::RePairReader<true> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(slp);

      auto report_compact_seq = [&compact_seq](const auto &_var) {
        compact_seq.emplace_back(_var);
      };

      re_pair_reader.Read(datafile, slp_wrapper);
    }

    cslp = grammar::CombinedSLP<>(slp);
    grammar::Chunks<> cslp_docs;
    grammar::Chunks<> cslp_occs;

    auto add_docs = [&cslp_docs, &cslp_occs](const auto &_slp, auto _var, auto ..._args) {
      auto set = _slp.Span(_var);
      std::map<std::size_t, std::size_t> counter;
      for (const auto &item : set) {
        ++counter[item];
      }

      std::vector<std::size_t> docs;
      docs.reserve(counter.size());
      std::vector<std::size_t> occs;
      occs.reserve(counter.size());

      for (const auto &item : counter) {
        docs.emplace_back(item.first);
        occs.emplace_back(item.second);
      }

      cslp_docs.Insert(docs.begin(), docs.end());
      cslp_occs.Insert(occs.begin(), occs.end());
    };

//    grammar::AddSet<decltype(cslp_docs)> add_set(cslp_docs);
    cslp.Compute(FLAGS_bs, add_docs, add_docs, grammar::MustBeSampled<decltype(cslp_docs)>(
        grammar::AreChildrenTooBig<decltype(cslp_docs)>(cslp_docs, FLAGS_sf)));

    sdsl::store_to_cache(cslp, KEY_GCDA_CSLP, config);

    sdsl::store_to_cache(cslp_docs, KEY_GCDA_CSLP_DOCS, config);
    grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
    sdsl::store_to_cache(cslp_docs_c, KEY_GCDA_CSLP_DOCS_C, config);

    sdsl::store_to_cache(cslp_occs, KEY_GCDA_CSLP_OCCS, config);
    grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_occs_c(cslp_occs, bit_compress, bit_compress);
    sdsl::store_to_cache(cslp_occs_c, KEY_GCDA_CSLP_OCCS_C, config);

    std::cout << "DONE" << std::endl;
  }

//  {
//    grammar::Chunks<> cslp_docs;
//    sdsl::load_from_cache(cslp_docs, KEY_GCDA_CSLP_DOCS, config);
//    grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
//    sdsl::store_to_cache(cslp_docs_c, KEY_GCDA_CSLP_DOCS_C, config);
//
//    grammar::Chunks<> cslp_occs;
//    sdsl::load_from_cache(cslp_occs, KEY_GCDA_CSLP_OCCS, config);
//    grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_occs_c(cslp_occs, bit_compress, bit_compress);
//    sdsl::store_to_cache(cslp_occs_c, KEY_GCDA_CSLP_OCCS_C, config);
//  }

  if (!sdsl::cache_file_exists(KEY_GCDA_LSLP, config)
      || !sdsl::cache_file_exists(KEY_GCDA_LSLP_BASIC, config)) {

    grammar::LightSLP<> lslp;
    if (!sdsl::load_from_cache(lslp, KEY_GCDA_LSLP, config)) {
      std::cout << "Construct LSLP" << std::endl;

      auto datafile = cache_file_name(KEY_DA_RAW, config);
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

      lslp.Compute(slp, compact_seq, cslp);

      sdsl::store_to_cache(lslp, KEY_GCDA_LSLP, config);

      std::cout << "DONE" << std::endl;
    }

    if (!sdsl::cache_file_exists(KEY_GCDA_LSLP_BASIC, config)) {
      std::cout << "Construct LSLP Basic" << std::endl;

      grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                        grammar::SampledSLP<>,
                        grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>
          lslp_basic(lslp, bit_compress, bit_compress, bit_compress, bit_compress);

      sdsl::store_to_cache(lslp_basic, KEY_GCDA_LSLP_BASIC, config);

      std::cout << "DONE" << std::endl;
    }
  }

  return 0;
}