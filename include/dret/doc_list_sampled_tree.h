//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <filesystem>

#include "sdsl/io.hpp"

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "construct_base.h"
#include "doc_list_index.h"
#include "index_base.h"

namespace dret {


template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
class DLSampledTreeScheme : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using Base = DocListIndexExtStorage<TStorage, TAlphabet>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  DLSampledTreeScheme(const TStorage& t_storage,
                      const TCountIdx& t_count_idx,
                      const TComputeCover& t_compute_cover,
                      const TGetDocs& t_get_docs,
                      const TGetDocSet& t_get_doc_set,
                      const TMergeSets& t_merge_sets)
      : Base(t_storage),
        count_idx_(t_count_idx),
        compute_cover_(t_compute_cover),
        get_docs_(t_get_docs),
        get_doc_set_(t_get_doc_set),
        merge_sets_(t_merge_sets) {}

  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = count_idx_.Count(t_pattern);

    auto cover = compute_cover_(sp, ep);

    const auto& range = cover.first;
    const auto& nodes = cover.second;

    std::vector<uint32_t> docs;
    docs.reserve(range.first - sp + ep - range.second);

    // TODO Use generic object to copy
    auto add_doc = [&docs](const auto& tt_d) {
      docs.emplace_back(tt_d);
    };

    if (nodes.empty()) {
      get_docs_(sp, ep, add_doc);
    } else {
      get_docs_(sp, range.first, add_doc);
      get_docs_(range.second, ep, add_doc);
    }

    sort(docs.begin(), docs.end());
    docs.erase(unique(docs.begin(), docs.end()), docs.end());

    if (!nodes.empty()) {
      merge_sets_(nodes.begin(), nodes.end(), get_doc_set_, docs);
    }

    for (const auto& doc : docs) {
      t_report(doc);
    }
  }

 protected:
  TCountIdx count_idx_;
  TComputeCover compute_cover_;
  TGetDocs get_docs_;
  TGetDocSet get_doc_set_;
  TMergeSets merge_sets_;
};

//~~~~~~~


template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
void construct(
    DLSampledTreeScheme<TStorage, TAlphabet, TCountIdx, TComputeCover, TGetDocs, TGetDocSet, TMergeSets>& t_index,
    Config& t_config) {
  if (!cache_file_exists(t_config.keys[conf::kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<TAlphabet::int_width>(t_config);
  }

  TCountIdx count_index(t_index.storage());
  construct(count_index, t_config.data_path, t_config);

  t_index.load(t_config);
}

//~~~~~~~


}  // namespace dret
