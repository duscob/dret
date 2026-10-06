//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/6/26.
//
// GCDA over a backend: GCDA's sampled tree and the document lists of its
// sampled nodes, with the partial leaves at the two ends of a range expanded by
// a document-array backend (plain DA, SA-Phi, RLCSA) instead of the grammar.
//
// The cover of a range reads only the sampled tree (Leaf, Position,
// IsFirstChild, Parent; see grammar::ComputeCoverFromBottom), so the grammar is
// needed to build the tree -- it is the sampled parse tree of the DA grammar --
// but not to query it. This index loads the cell's da_sampled_tree and node
// lists components and no grammar file. Against PDL over the same backend it
// isolates the sampling structure: both expand at most ~2b positions at the
// ends of a range with the same backend, and differ only in where the stored
// lists sit (a sample of the DA grammar's parse tree vs of the suffix tree).
//

#pragma once

#include <functional>
#include <utility>
#include <vector>

#include <grammar/sampled_slp.h>
#include <grammar/slp_helper.h>

#include "sr-index/r_index.h"

#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_sampled_tree_base.h"
#include "dret/index_base.h"
#include "dret/pdl/get_docs.h"
#include "dret/slp/slp_components.h"

namespace dret::gcda {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TGetDocs = pdl::PDLGetDocsSAPhi_R<TStorage, TAlphabet::int_width>,
          typename TSets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                             true,
                                             grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TSampledTree = grammar::SampledSLP<>,
          typename TMergeSets = MergeSetsBinaryTreeFunctor>
class DocListIdxGCDABackend : public DLSampledTreeScheme<TMergeSets, typename TAlphabet::string_type>,
                              public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  using SchemeBase = DLSampledTreeScheme<TMergeSets, typename TAlphabet::string_type>;
  using StorageBase = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using typename SchemeBase::size_type;
  using GetDocs = TGetDocs;
  using Sets = TSets;

  explicit DocListIdxGCDABackend(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_storage),
        get_docs_(typename TGetDocs::Inner{t_storage}),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor) {}

  DocListIdxGCDABackend() = default;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    std::size_t bytes = count_idx_.serialize(out, child, "count_idx");
    bytes += get_docs_.inner().serialize(out, child, "get_docs");
    if (tree_) bytes += tree_->serialize(out, child, "sampled_tree");
    if (sets_) bytes += sets_->serialize(out, child, "sets");
    return bytes;
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    for (const auto& f : get_docs_.inner().GetSizeReport()) append(r, "get_docs_" + f.name, f.bytes);
    if (tree_) append(r, "sampled_tree", sdsl::size_in_bytes(*tree_));
    if (sets_) append(r, "sets", sdsl::size_in_bytes(*sets_));
    return r;
  }

  const TCountIdx& count_idx() const { return count_idx_; }
  TGetDocs& get_docs() { return get_docs_; }
  uint32_t block_size() const { return block_size_; }
  float storing_factor() const { return storing_factor_; }

 protected:
  std::pair<std::size_t, std::size_t> count(const typename SchemeBase::TPattern& t_pattern) const override {
    return count_idx_.Count(t_pattern);
  }

  std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> computeCover(
      std::size_t t_sp,
      std::size_t t_ep) const override {
    std::vector<std::size_t> nodes;
    auto report = [&nodes](const auto& _value) { nodes.emplace_back(_value); };
    auto range = grammar::ComputeCoverFromBottom(*tree_, t_sp, t_ep, report);
    return {std::move(range), std::move(nodes)};
  }

  void getDocs(std::size_t t_sp,
               std::size_t t_ep,
               const std::function<void(typename SchemeBase::TDocId)>& t_report) const override {
    auto report = [&t_report](std::size_t d) { t_report(d); };
    get_docs_.getDocs(t_sp, t_ep, report);
  }

  std::vector<typename SchemeBase::TDocId> getDocSet(std::size_t t_i) const override {
    auto v = (*sets_)[t_i];
    return {v.begin(), v.end()};
  }

  void loadInner(typename StorageBase::TSource& t_source, const JSON& t_keys) override {
    std::visit(
        [this](auto&& tt_source) {
          count_idx_.load(tt_source.get());
          get_docs_.inner().load(tt_source.get());
        },
        t_source);
    const SampledTreeCell cell{block_size_, storing_factor_};
    tree_ = this->template loadItemPtr<TSampledTree>(CellKey(t_keys, conf::kDaSampledTree, cell), t_source, true);
    sets_ = this->template loadItemPtr<TSets>(NodeDocListsKey<TSets>(t_keys, cell), t_source, true);
  }

  TCountIdx count_idx_;
  TGetDocs get_docs_;
  const TSampledTree* tree_ = nullptr;
  const TSets* sets_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

// The sampled tree and the node lists come from the cell's GCDA build (which
// reads the DA grammar); the backend builds its own cache.
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TGetDocs, typename TSets,
          typename TSampledTree, typename TMergeSets>
void construct(DocListIdxGCDABackend<TStorage, TAlphabet, TCountIdx, TGetDocs, TSets, TSampledTree, TMergeSets>& t_index,
               Config& t_config) {
  using namespace conf;
  if (!cache_file_exists(t_config.keys[kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<TAlphabet::int_width>(t_config);
  }
  if (const auto key = t_config.keys[kSA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    sdsl::construct_sa<TAlphabet::int_width>(t_config);
  }
  if (const auto key = t_config.keys[kDocEnds].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocEnd<TAlphabet::int_width, sdsl::sd_vector<>>(t_config);
  }
  if (const auto key = t_config.keys[kDA].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray(t_config);
  }

  const SampledTreeCell cell{t_index.block_size(), t_index.storing_factor()};
  if (!sdsl::cache_file_exists<TSampledTree>(CellKey(t_config.keys, kDaSampledTree, cell), t_config)) {
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
    LoadOrBuildSampledTree(t_config, filepath_da, cell.block_size, cell.storing_factor);
  }
  if (const auto key = NodeDocListsKey<TSets>(t_config.keys, cell); !sdsl::cache_file_exists<TSets>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    TSets sets;
    construct(sets, t_config, cell.block_size, cell.storing_factor);
  }

  construct(t_index.get_docs(), t_config);
  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  t_index.load(t_config);
}

}  // namespace dret::gcda
