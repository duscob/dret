//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <algorithm>
#include <functional>
#include <vector>

#include "dret/doc_list/doc_list_base.h"

#ifdef DRET_DOC_LIST_PROFILE
#include <chrono>

#include "dret/doc_list/search_profile.h"
#endif

namespace dret {


template <typename TMergeSets, typename TSequence = Alphabet<>::string_type>
class DLSampledTreeScheme : public DocListIndex<TSequence> {
 public:
  using Base = DocListIndex<TSequence>;
  using typename Base::TDocId;
  using typename Base::TPattern;
  using size_type = std::size_t;

  explicit DLSampledTreeScheme(const TMergeSets& t_merge_sets) : merge_sets_(t_merge_sets) {}

  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
#ifdef DRET_DOC_LIST_PROFILE
    using prof_clk = std::chrono::steady_clock;
    const auto prof_t0 = prof_clk::now();
#endif
    auto [sp, ep] = count(t_pattern);
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t1 = prof_clk::now();
#endif

    std::vector<std::pair<std::size_t, std::size_t>> raw_ranges;
    std::vector<std::size_t> nodes;
    computeCoverFull(sp, ep, raw_ranges, nodes);
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t2 = prof_clk::now();
#endif

    std::size_t total_raw = 0;
    for (const auto& r : raw_ranges) total_raw += r.second - r.first;
    std::vector<TDocId> docs;
    docs.reserve(total_raw);

    auto add_doc = [&docs](const auto& tt_d) {
      docs.emplace_back(tt_d);
    };

    for (const auto& [b, e] : raw_ranges) {
      getDocs(b, e, add_doc);
    }

    sort(docs.begin(), docs.end());
    docs.erase(unique(docs.begin(), docs.end()), docs.end());
#ifdef DRET_DOC_LIST_PROFILE
    const auto prof_t3 = prof_clk::now();
#endif

    if (!nodes.empty()) {
      auto get_doc_set = [this](std::size_t t_i) {
        return this->getDocSet(t_i);
      };
      merge_sets_(nodes.begin(), nodes.end(), get_doc_set, docs);
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
      pr.n_raw_positions += total_raw;
      pr.n_docs += docs.size();
    }
#endif

    for (const auto& doc : docs) {
      t_report(doc);
    }
  }

 protected:
  virtual std::pair<std::size_t, std::size_t> count(const TPattern& t_pattern) const = 0;

  virtual std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> computeCover(
      std::size_t t_sp,
      std::size_t t_ep) const = 0;

  // Multi-range cover. Default delegates to the single-range computeCover
  // and produces at most a prefix and a suffix raw range — semantically
  // equivalent to the previous Search shape, so GCDA/DGCDA keep working
  // unchanged. PDL indexes override this to return interleaved raw gaps
  // when their storage policy leaves leaves unselected (see
  // PDLTreeCore::computeCoverFull, Task 12).
  virtual void computeCoverFull(std::size_t t_sp,
                                std::size_t t_ep,
                                std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                                std::vector<std::size_t>& t_nodes) const {
    auto [range, nds] = computeCover(t_sp, t_ep);
    if (nds.empty()) {
      t_raw_ranges.emplace_back(t_sp, t_ep);
    } else {
      if (t_sp < range.first) t_raw_ranges.emplace_back(t_sp, range.first);
      if (range.second < t_ep) t_raw_ranges.emplace_back(range.second, t_ep);
    }
    t_nodes = std::move(nds);
  }

  virtual void getDocs(std::size_t t_sp, std::size_t t_ep, const std::function<void(TDocId)>& t_report) const = 0;

  virtual std::vector<TDocId> getDocSet(std::size_t t_i) const = 0;

  TMergeSets merge_sets_;
};

}  // namespace dret
