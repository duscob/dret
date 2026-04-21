//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <algorithm>
#include <functional>
#include <vector>

#include "doc_list_index.h"

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
    auto [sp, ep] = count(t_pattern);

    auto cover = computeCover(sp, ep);

    const auto& range = cover.first;
    const auto& nodes = cover.second;

    std::vector<uint32_t> docs;
    docs.reserve(range.first - sp + ep - range.second);

    auto add_doc = [&docs](const auto& tt_d) {
      docs.emplace_back(tt_d);
    };

    if (nodes.empty()) {
      getDocs(sp, ep, add_doc);
    } else {
      getDocs(sp, range.first, add_doc);
      getDocs(range.second, ep, add_doc);
    }

    sort(docs.begin(), docs.end());
    docs.erase(unique(docs.begin(), docs.end()), docs.end());

    if (!nodes.empty()) {
      auto get_doc_set = [this](std::size_t t_i) { return this->getDocSet(t_i); };
      merge_sets_(nodes.begin(), nodes.end(), get_doc_set, docs);
    }

    for (const auto& doc : docs) {
      t_report(doc);
    }
  }

 protected:
  virtual std::pair<std::size_t, std::size_t> count(const TPattern& t_pattern) const = 0;

  virtual std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>>
      computeCover(std::size_t t_sp, std::size_t t_ep) const = 0;

  virtual void getDocs(std::size_t t_sp,
                       std::size_t t_ep,
                       const std::function<void(std::size_t)>& t_report) const = 0;

  virtual std::vector<uint32_t> getDocSet(std::size_t t_i) const = 0;

  TMergeSets merge_sets_;
};

//~~~~~~~


}  // namespace dret
