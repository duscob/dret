//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//
// The core of the sampled-tree indexes (GCDA, GCDA over a backend, PDL): cover the suffix-array range
// with the sampled nodes inside it, read the stored lists of those nodes, and expand the positions no
// node covers through a document-array backend (backend/backend.h), every document into a sink
// (dret/sink.h). The range itself comes from the index (doc_list_index.h).
//
// The indexes differ only in the sampler, which owns the sampled nodes and their lists:
//
//   GrammarTree  the sampled parse tree of the DA grammar, read from the backend's grammar, which also
//                expands the uncovered ends (GCDA; doc_list_gcda.h)
//   SampledTree  the same tree stored alone, with any backend for the ends (GCDA over a backend;
//                doc_list_gcda_backend.h)
//   PDLTree      a sample of the suffix tree (PDL; doc_list_pdl.h)
//
// A sampler is an IndexBaseWithExternalStorage (it loads its structures from the shared storage),
// constructed from (storage, block size, storing factor, ...), with
//
//   computeCoverFull(backend, sp, ep, raw, nodes)  the uncovered ranges and the covering nodes
//   forEachInNode(node, f)                         f(d) for every document of the node's list
//   kBackendSharesCell                             the backend is built from the sampler's block size
//                                                  and storing factor (GCDA's grammar), else from the
//                                                  storage alone, unless the caller sets it up
//   kBackendName, kBackendPrefix                   the backend's name in the serialization and the
//                                                  prefix of its fields in the size report
//   serialize, GetSizeReport, block_size, storing_factor
//
// and a construct(sampler, config) that builds its structures.
//
// Serialized as the backend, then the sampler: the layout each of the three indexes had on its own.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <sdsl/io.hpp>

#include "dret/backend/backend.h"
#include "dret/collection.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

#ifdef DRET_DOC_LIST_PROFILE
#include <chrono>

#include "dret/doc_list/search_profile.h"
#endif

namespace dret {

// The uncovered ranges of a single-range cover: all of [sp, ep) when no node covers any of it,
// otherwise the parts before and after the span [range.first, range.second) the nodes cover.
inline void CoverToRanges(std::size_t t_sp,
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

template <typename TStorage, uint8_t t_width, typename TBackend, typename TSampler>
class SampledTreeCore
    : public FieldsComponent<SampledTreeCore<TStorage, t_width, TBackend, TSampler>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<SampledTreeCore, TStorage, t_width>;
  using typename Base::size_type;
  using Backend = TBackend;
  using Sampler = TSampler;

  // The query profile (search_profile.h) also times the index's count.
  static constexpr bool kProfiled = true;

  // The backend from the storage (and, for GCDA's grammar, the sampler's cell); t_sampler_args after
  // the storing factor go to the sampler (PDL's storage policy).
  template <typename... TArgs>
  explicit SampledTreeCore(const TStorage& t_storage,
                           uint32_t t_block_size = 512,
                           float t_storing_factor = 4,
                           TArgs&&... t_sampler_args)
      : Base(t_storage),
        backend_(MakeBackend(t_storage, t_block_size, t_storing_factor)),
        sampler_(t_storage, t_block_size, t_storing_factor, std::forward<TArgs>(t_sampler_args)...) {}

  // A backend set up by the caller, with its own parameters (e.g. a grammar's cell for PDL).
  template <typename... TArgs>
  SampledTreeCore(const TStorage& t_storage,
                  const TBackend& t_backend,
                  uint32_t t_block_size = 512,
                  float t_storing_factor = 4,
                  TArgs&&... t_sampler_args)
      : Base(t_storage),
        backend_(t_backend),
        sampler_(t_storage, t_block_size, t_storing_factor, std::forward<TArgs>(t_sampler_args)...) {}

  SampledTreeCore() = default;

  template <typename TSink>
  void Query(std::size_t t_sp, std::size_t t_ep, std::size_t /*t_m*/, TSink& t_sink) const {
#ifdef DRET_DOC_LIST_PROFILE
    using prof_clk = std::chrono::steady_clock;
    const auto prof_t1 = prof_clk::now();
#endif
    std::vector<std::pair<std::size_t, std::size_t>> raw_ranges;
    std::vector<std::size_t> nodes;
    sampler_.computeCoverFull(backend_, t_sp, t_ep, raw_ranges, nodes);
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t2 = prof_clk::now();
#endif

    // A reference to the sink, not the sink: some expansions take their report as const&.
    auto position = [&t_sink](auto d) { t_sink.Position(static_cast<std::size_t>(d)); };
    for (const auto& [b, e] : raw_ranges) {
      if (t_sink.Done()) return;
      backend_(b, e, position);
    }
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t3 = prof_clk::now();
#endif

    auto listed = [&t_sink](std::size_t d) { t_sink.Listed(d); };
    for (auto i : nodes) {
      if (t_sink.Done()) return;
      sampler_.forEachInNode(i, listed);
    }
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t4 = prof_clk::now();
    {
      using ns = std::chrono::nanoseconds;
      auto& pr = dret::search_profile();
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

  // A backend that finds the range itself (the RLCSA, as in the original PDL) finds it for the index,
  // which then holds no r-index.
  template <typename TPattern>
    requires ::dret::backend::CountsPattern<TBackend, TPattern>
  std::pair<std::size_t, std::size_t> Count(const TPattern& t_pattern) const {
    return backend_.Count(t_pattern);
  }

  // Not stored: the listing bitmap grows as needed (distinct_docs.h).
  static constexpr std::size_t n_doc() { return 0; }

  const TBackend& backend() const { return backend_; }

  const TSampler& sampler() const { return sampler_; }

  uint32_t block_size() const { return sampler_.block_size(); }

  float storing_factor() const { return sampler_.storing_factor(); }

 protected:
  friend Base;

  static TBackend MakeBackend(const TStorage& t_storage, uint32_t t_block_size, float t_storing_factor) {
    if constexpr (TSampler::kBackendSharesCell)
      return TBackend(t_storage, t_block_size, t_storing_factor);
    else
      return TBackend(t_storage);
  }

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    return std::tuple{Part(TSampler::kBackendName, t_self.backend_, TSampler::kBackendPrefix),
                      Part("sampler", t_self.sampler_)};
  }

  TBackend backend_;
  TSampler sampler_;
};

// The DA, then the sampler's structures and the backend's: the backend first when it is the grammar the
// sampler reads its tree from (GCDA), as each index built them on its own. Copies, as neither need be
// default-constructible (not when TStorage is std::reference_wrapper).
template <typename TStorage, uint8_t t_width, typename TBackend, typename TSampler>
void construct(SampledTreeCore<TStorage, t_width, TBackend, TSampler>& t_core, Config& t_config) {
  EnsureDA<t_width>(t_config);

  auto backend = t_core.backend();
  auto sampler = t_core.sampler();
  if constexpr (TSampler::kBackendSharesCell) {
    construct(backend, t_config);
    construct(sampler, t_config);
  } else {
    construct(sampler, t_config);
    construct(backend, t_config);
  }
}

}  // namespace dret
