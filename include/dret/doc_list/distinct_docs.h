//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/19/26.
//
// The deduplication shared by the brute-force indexes, which read the document of every position of
// the suffix-array range from a backend (SA-Phi, the sr-index, the plain DA, the RLCSA, a GCDA
// grammar) and report each document once.
//
// A bitmap over the document identifiers, sized to the number of documents and cleared per query:
// a document is reported the first time it is seen, so a query costs O(occ + n_doc / 64) and reports
// in discovery order (set semantics; see DocListResult in bm_doc_list.cpp). Sorting the occurrences
// instead costs O(occ log occ), which is what the many occurrences per document of the repetitive
// collections make expensive. Identifiers past n_doc grow the bitmap rather than being dropped.
//
// The number of documents is cached at conf::kRmqNDoc, the entry the RMQ cores read.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>

#include "dret/config.h"

namespace dret {

class DistinctDocs {
 public:
  explicit DistinctDocs(std::size_t t_n_doc) : words_((t_n_doc >> 6) + 1, 0) {}

  // True the first time t_doc is inserted.
  bool Insert(std::size_t t_doc) {
    const std::size_t w = t_doc >> 6;
    if (w >= words_.size()) words_.resize(w + 1, 0);
    const uint64_t bit = uint64_t{1} << (t_doc & 63);
    if (words_[w] & bit) return false;
    words_[w] |= bit;
    return true;
  }

 private:
  std::vector<uint64_t> words_;
};

// Caches the number of documents (the kInternalDocDelim bytes of the text) at conf::kRmqNDoc, unless
// it is there already.
inline void ConstructNDoc(Config& t_config) {
  const auto key = t_config.keys[conf::kRmqNDoc].get<std::string>();
  if (sdsl::cache_file_exists<sdsl::int_vector<64>>(key, t_config)) return;

  sdsl::int_vector_buffer<8> text_buf(sdsl::cache_file_name(t_config.keys[conf::kText].get<std::string>(), t_config));
  std::size_t n_doc = 0;
  for (std::size_t i = 0; i < text_buf.size(); ++i) {
    if (text_buf[i] == 1) ++n_doc;
  }
  sdsl::int_vector<64> n_doc_vec(1, n_doc);
  sdsl::store_to_cache(n_doc_vec, key, t_config, true);
}

}  // namespace dret
