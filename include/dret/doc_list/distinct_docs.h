//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/19/26.
//
// The document bitmap behind ListSink (doc_list/list_sink.h), which every listing index reports through: each
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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace dret {

// The words and their count are plain members, not a std::vector: Insert's bound check then reads one
// scalar. With a vector, the check is (end - begin) / 8, and GCC with AVX-512 (-march=native on shannon)
// packed the sink's pointers into one vector register and computed the size with vector instructions on
// every position, 10% on Brute-DA over Revision.
class DistinctDocs {
 public:
  explicit DistinctDocs(std::size_t t_n_doc)
      : n_words_((t_n_doc >> 6) + 1), words_(std::make_unique<uint64_t[]>(n_words_)) {}

  // True the first time t_doc is inserted. Forced inline into every report (GCC 14 on shannon left it
  // out of line in the ILCP-family cores, and VisitRun with it: 6% on CILCP over the DA); growing the
  // bitmap, which happens only past the size given, is out of line.
  [[gnu::always_inline]] bool Insert(std::size_t t_doc) {
    const std::size_t w = t_doc >> 6;
    if (w >= n_words_) [[unlikely]] Grow(w);
    const uint64_t bit = uint64_t{1} << (t_doc & 63);
    uint64_t& word = words_[w];
    if (word & bit) return false;
    word |= bit;
    return true;
  }

  // Whether t_doc was inserted. Forced inline like Insert: the RMQ cores' marking stop calls it.
  [[gnu::always_inline]] bool Contains(std::size_t t_doc) const {
    const std::size_t w = t_doc >> 6;
    return w < n_words_ && (words_[w] & (uint64_t{1} << (t_doc & 63)));
  }

 private:
  [[gnu::noinline]] void Grow(std::size_t t_word) {
    auto words = std::make_unique<uint64_t[]>(t_word + 1);
    std::copy(words_.get(), words_.get() + n_words_, words.get());
    words_ = std::move(words);
    n_words_ = t_word + 1;
  }

  std::size_t n_words_;
  std::unique_ptr<uint64_t[]> words_;
};

}  // namespace dret
