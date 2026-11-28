//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// Result sinks: what a search hands its documents to. A core finds the documents of a suffix-array
// range -- by expanding positions through a backend, by reading the stored list of a sampled node,
// by visiting the runs of an RMQ partition -- and a sink decides what the query returns. The cores
// are shared by every task; the sink is the task.
//
//   Position(d)  d is the document of one position of the range
//   Listed(d)    d occurs in the range (a stored list: how often is not known)
//   Done()       the search may stop (a sink that has its answer, e.g. top-k with a threshold)
//
// A sink is also callable as a position report, operator()(d) == Position(d), so it can be handed
// to a backend's range expansion directly.
//
// ListSink is document listing: each document reported once, the first time it is seen. FreqSink is
// document frequencies: how many positions of the range each document has, from the positions or
// from stored lists that carry counts,
//
//   Weighted(d, f)  d occurs f times in a part of the range (a stored list with counts)
//   Finish()        the search is over: report each document with its frequency
//
// A core serves frequencies when it hands the sink every position of the range once, or the counts
// of the parts it does not expand (kFrequencies, doc_list_index.h). A list without counts cannot
// serve it: FreqSink::Listed does not compile.
//

#pragma once

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "dret/doc_list/distinct_docs.h"

namespace dret {

template <typename TSink>
concept DocSink = requires(TSink& t_sink, std::size_t t_doc) {
  t_sink.Position(t_doc);
  t_sink.Listed(t_doc);
  { t_sink.Done() } -> std::convertible_to<bool>;
};

// Reports t_doc, found at suffix-array position t_pos, to t_sink: with the position when the sink takes
// it (PositionAt: the extreme positions the frequency cores collect, doc_freq/rmq_core.h), else as a
// position. The RMQ cores report through it.
template <typename TSink>
[[gnu::always_inline]] inline void ReportPosition(TSink& t_sink, std::size_t t_doc, std::size_t t_pos) {
  if constexpr (requires { t_sink.PositionAt(t_doc, t_pos); })
    t_sink.PositionAt(t_doc, t_pos);
  else
    t_sink.Position(t_doc);
}

// Document listing: reports each document once, to t_report, the first time it is seen, from a
// position or a list. Contains(d) tells whether d was reported, for the marking stop of the RMQ cores.
// t_n_doc sizes the bitmap; identifiers past it grow it (distinct_docs.h), so 0 is valid when the
// number of documents is not at hand.
template <typename TReport>
class ListSink {
 public:
  ListSink(std::size_t t_n_doc, const TReport& t_report) : seen_(t_n_doc), report_(t_report) {}

  // Inlined into the cores' inner loops (an RMQ run's fan-out reports every position): a call per
  // position measured 6% on ILCP over the DA (revision-small).
  [[gnu::always_inline]] void Position(std::size_t t_doc) {
    if (seen_.Insert(t_doc)) {
#ifdef DRET_DOC_LIST_PROFILE
      ++reported_;
#endif
      report_(t_doc);
    }
  }

  [[gnu::always_inline]] void Listed(std::size_t t_doc) { Position(t_doc); }

  [[gnu::always_inline]] void operator()(std::size_t t_doc) { Position(t_doc); }

  static constexpr bool Done() { return false; }

  bool Contains(std::size_t t_doc) const { return seen_.Contains(t_doc); }

#ifdef DRET_DOC_LIST_PROFILE
  // Documents reported, for the query profile (dret/doc_list/search_profile.h).
  std::size_t Reported() const { return reported_; }
#endif

 private:
  DistinctDocs seen_;
  const TReport& report_;
#ifdef DRET_DOC_LIST_PROFILE
  std::size_t reported_ = 0;
#endif
};

// Document frequencies: adds up, per document, the positions and the counts the core hands it, and
// reports each document once with its total, in Finish(), in the order the documents were first seen.
// t_n_doc sizes the counters; identifiers past it grow them, so 0 is valid as for ListSink.
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

}  // namespace dret
