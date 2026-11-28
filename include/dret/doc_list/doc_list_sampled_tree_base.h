//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//
// The query of the sampled-tree indexes (GCDA, GCDA over a backend, PDL): find the suffix-array range
// of the pattern, cover it with the sampled nodes inside it, read the stored lists of those nodes and
// expand the positions no node covers through the index's backend, every document into a sink
// (dret/sink.h).
//
// The index supplies the steps as members, found at compile time (CRTP), so a position reaches the
// sink without a virtual call or a std::function:
//
//   count(pattern)                          -> [sp, ep), the suffix-array range
//   computeCoverFull(sp, ep, raw, nodes)    the uncovered ranges and the covering nodes
//   getDocs(b, e, sink)                     sink.Position(d) for every position of [b, e)
//   forEachInNode(node, sink)               sink.Listed(d) for every document of the node's list
//
// CoverToRanges turns a single-range cover (the nodes and the span they cover) into the uncovered
// ranges, for the indexes whose cover is one contiguous span (GCDA).
//

#pragma once

#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "dret/doc_list/doc_list_base.h"
#include "dret/sink.h"

#ifdef DRET_DOC_LIST_PROFILE
#include <chrono>

#include "dret/doc_list/search_profile.h"
#endif

namespace dret {

template <typename TDerived, typename TMergeSets, typename TSequence = Alphabet<>::string_type>
class DLSampledTreeScheme : public DocListIndex<TSequence> {
 public:
  using Base = DocListIndex<TSequence>;
  using typename Base::TDocId;
  using typename Base::TPattern;
  using size_type = std::size_t;

  explicit DLSampledTreeScheme(const TMergeSets& t_merge_sets) : merge_sets_(t_merge_sets) {}

  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    // Each document is reported the first time it is seen, from the uncovered ranges or from a stored
    // list. The number of documents is not at hand here, so the bitmap grows as needed.
    ListSink<std::function<void(TDocId)>> sink(0, t_report);
    Query(t_pattern, sink);
  }

  // The query into any sink.
  template <typename TSink>
  void Query(const TPattern& t_pattern, TSink& t_sink) const {
    const auto& self = static_cast<const TDerived&>(*this);
#ifdef DRET_DOC_LIST_PROFILE
    using prof_clk = std::chrono::steady_clock;
    const auto prof_t0 = prof_clk::now();
#endif
    auto [sp, ep] = self.count(t_pattern);
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t1 = prof_clk::now();
#endif

    std::vector<std::pair<std::size_t, std::size_t>> raw_ranges;
    std::vector<std::size_t> nodes;
    self.computeCoverFull(sp, ep, raw_ranges, nodes);
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t2 = prof_clk::now();
#endif

    for (const auto& [b, e] : raw_ranges) {
      if (t_sink.Done()) return;
      self.getDocs(b, e, t_sink);
    }
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t3 = prof_clk::now();
#endif

    for (auto i : nodes) {
      if (t_sink.Done()) return;
      self.forEachInNode(i, t_sink);
    }
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t4 = prof_clk::now();
    {
      using ns = std::chrono::nanoseconds;
      auto& pr = dret::search_profile();
      pr.ns_count += std::chrono::duration_cast<ns>(prof_t1 - prof_t0).count();
      pr.ns_cover += std::chrono::duration_cast<ns>(prof_t2 - prof_t1).count();
      pr.ns_expand += std::chrono::duration_cast<ns>(prof_t3 - prof_t2).count();
      pr.ns_combine += std::chrono::duration_cast<ns>(prof_t4 - prof_t3).count();
      ++pr.n_queries;
      pr.n_nodes += nodes.size();
      for (const auto& r : raw_ranges) pr.n_raw_positions += r.second - r.first;
      if constexpr (requires { t_sink.Reported(); }) pr.n_docs += t_sink.Reported();
    }
#endif
  }

 protected:
  // The uncovered ranges of a single-range cover: all of [sp, ep) when no node covers any of it,
  // otherwise the parts before and after the span [range.first, range.second) the nodes cover.
  static void CoverToRanges(std::size_t t_sp,
                            std::size_t t_ep,
                            const std::pair<std::size_t, std::size_t>& t_range,
                            bool t_has_nodes,
                            std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges) {
    if (!t_has_nodes) {
      t_raw_ranges.emplace_back(t_sp, t_ep);
      return;
    }
    if (t_sp < t_range.first) t_raw_ranges.emplace_back(t_sp, t_range.first);
    if (t_range.second < t_ep) t_raw_ranges.emplace_back(t_range.second, t_ep);
  }

  // Not used by Query since the bitmap deduplication (distinct_docs.h); kept for the
  // MergeSetsBinaryTreeFunctor of the index types and its tests.
  TMergeSets merge_sets_;
};

}  // namespace dret
