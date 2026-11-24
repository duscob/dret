//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/19/26.
//
// The document bitmap behind ListSink (dret/sink.h), which every listing index reports through: each
// document reported once, whether it comes from a position of the range, a stored list or an RMQ run.
//
// A bitmap over the document identifiers, sized to the number of documents and cleared per query:
// a document is reported the first time it is seen, so a query costs O(occ + n_doc / 64) and reports
// in discovery order (set semantics; see DocListResult in bm_doc_list.cpp). Sorting the occurrences
// instead costs O(occ log occ), which is what the many occurrences per document of the repetitive
// collections make expensive. Identifiers past n_doc grow the bitmap rather than being dropped.
//
// The number of documents comes from dret/collection.h (EnsureNDoc).
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

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

  // Whether t_doc was inserted.
  bool Contains(std::size_t t_doc) const {
    const std::size_t w = t_doc >> 6;
    return w < words_.size() && (words_[w] & (uint64_t{1} << (t_doc & 63)));
  }

 private:
  std::vector<uint64_t> words_;
};

}  // namespace dret
