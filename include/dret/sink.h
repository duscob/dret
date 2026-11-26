//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// Result sinks: what a search hands its documents to. A core finds the documents of a suffix-array
// range -- by expanding positions through a backend, by reading the stored list of a sampled node,
// by visiting the runs of an RMQ partition -- and a sink decides what the query returns. The cores
// are shared by every task; the sink is the task. The contract:
//
//   Position(d)       d is the document of one position of the range
//   Listed(d)         d occurs in the range (a stored list: how often is not known)
//   Done()            the search may stop (a sink that has its answer, e.g. top-k with a threshold)
//
// and, for the sinks that want more:
//
//   PositionAt(d, i)  d at suffix-array position i (the cores that know the position hand it over
//                     through ReportPosition below)
//   Weighted(d, f)    d occurs f times in a part of the range (a stored list with counts)
//
// A sink is also callable as a position report, operator()(d) == Position(d), so it can be handed
// to a backend's range expansion directly.
//
// The sinks themselves belong to their tasks: ListSink to listing (doc_list/list_sink.h), FreqSink to
// frequencies (doc_freq/frequencies.h).
//

#pragma once

#include <concepts>
#include <cstddef>

namespace dret {

template <typename TSink>
concept DocSink = requires(TSink& t_sink, std::size_t t_doc) {
  t_sink.Position(t_doc);
  t_sink.Listed(t_doc);
  { t_sink.Done() } -> std::convertible_to<bool>;
};

// Reports t_doc, found at suffix-array position t_pos, to t_sink: with the position when the sink takes
// it (PositionAt), else as a position. The RMQ cores report through it.
template <typename TSink>
[[gnu::always_inline]] inline void ReportPosition(TSink& t_sink, std::size_t t_doc, std::size_t t_pos) {
  if constexpr (requires { t_sink.PositionAt(t_doc, t_pos); })
    t_sink.PositionAt(t_doc, t_pos);
  else
    t_sink.Position(t_doc);
}

}  // namespace dret
