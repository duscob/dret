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

  DLSampledTreeScheme(const TStorage& t_storage) : Base(t_storage), count_idx_(t_storage) {}

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

  void load(Config t_config) override {
    count_idx_.load(t_config);
    compute_cover_.load(t_config);
    get_docs_.load(t_config);
    get_doc_set_.load(t_config);
    merge_sets_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    count_idx_.load(in);
    compute_cover_.load(in);
    get_docs_.load(in);
    get_doc_set_.load(in);
    merge_sets_.load(in);
  }

  using Base::load;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return count_idx_.serialize(out, child, "count_idx") + compute_cover_.serialize(out, child, "compute_cover")
           + get_docs_.serialize(out, child, "get_docs") + get_doc_set_.serialize(out, child, "get_doc_set")
           + merge_sets_.serialize(out, child, "merge_sets");
  }

 protected:
  TCountIdx count_idx_;
  TComputeCover compute_cover_;
  const TGetDocs get_docs_;
  const TGetDocSet get_doc_set_;
  const TMergeSets merge_sets_;
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
