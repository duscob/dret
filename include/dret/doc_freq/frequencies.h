//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Document frequencies over the listing indexes: the frequency sink (the contract is dret/sink.h),
// which cores serve it, and the query.
//
// A core serves frequencies when it hands the sink every position of the range once (brute force), or
// the counts of the parts it does not expand (GCDA with counts, doc_freq/gcda_counts.h), or the
// documents with their frequencies (the cores with an oracle, doc_freq/rmq_core.h and
// doc_freq/gcda_occs_core.h). The cores are listing's (doc_list/); which of them serve frequencies is
// said here and next to the frequency components, by specializing ComputesFrequencies, so the listing
// code carries nothing of it.
//

#pragma once

#include <cstddef>
#include <type_traits>
#include <vector>

#include "dret/doc_list/brute_core.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/sink.h"

namespace dret::freq {

// Document frequencies: adds up, per document, the positions and the counts the core hands it, and
// reports each document once with its total, in Finish(), in the order the documents were first seen.
// t_n_doc sizes the counters; identifiers past it grow them, so 0 is valid.
template <typename TReport>
class FreqSink {
 public:
  FreqSink(std::size_t t_n_doc, const TReport& t_report) : counts_(t_n_doc, 0), report_(t_report) {}

  [[gnu::always_inline]] void Position(std::size_t t_doc) { Weighted(t_doc, 1); }

  [[gnu::always_inline]] void operator()(std::size_t t_doc) { Position(t_doc); }

  [[gnu::always_inline]] void Weighted(std::size_t t_doc, std::size_t t_freq) {
    if (t_doc >= counts_.size()) counts_.resize(t_doc + 1, 0);
    if (counts_[t_doc] == 0) docs_.push_back(t_doc);
    counts_[t_doc] += t_freq;
  }

  // A stored list without counts says a document occurs, not how often.
  template <typename T = void>
  void Listed(std::size_t) {
    static_assert(!std::is_same_v<T, T>, "FreqSink: this core reads stored lists without counts, so it cannot compute "
                                  "document frequencies");
  }

  static constexpr bool Done() { return false; }

  void Finish() {
    for (auto d : docs_) {
      report_(d, counts_[d]);
      counts_[d] = 0;
    }
    docs_.clear();
  }

 private:
  std::vector<std::size_t> counts_;
  std::vector<std::size_t> docs_;  // the documents with a count, in the order first seen
  const TReport& report_;
};

// Whether a core's Query serves FreqSink.
template <typename TCore>
struct ComputesFrequencies : std::false_type {};

// Brute force hands every position over once.
template <typename TStorage, uint8_t t_width, typename TBackend>
struct ComputesFrequencies<BruteCore<TStorage, t_width, TBackend>> : std::true_type {};

// The locate cores hand every occurrence over once.
template <typename TStorage, typename TAlphabet, typename TLocateIdx, typename TGetDoc>
struct ComputesFrequencies<LocateCore<TStorage, TAlphabet, TLocateIdx, TGetDoc>> : std::true_type {};

template <typename TIndex>
concept FrequencyIndex = ComputesFrequencies<typename TIndex::Core>::value;

// The documents of t_pattern with their frequencies: t_report(d, f) once per document, in no given
// order.
template <typename TIndex, typename TPattern, typename TReport>
  requires FrequencyIndex<TIndex>
void Frequencies(const TIndex& t_index, const TPattern& t_pattern, const TReport& t_report) {
  FreqSink<TReport> sink(t_index.core().n_doc(), t_report);
  t_index.Query(t_pattern, sink);
  sink.Finish();
}

}  // namespace dret::freq
