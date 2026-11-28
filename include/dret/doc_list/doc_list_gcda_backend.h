//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/6/26.
//
// GCDA over a backend: the sampled-tree core (sampled_tree_core.h) with GCDA's
// sampled tree and the document lists of its sampled nodes (SampledTree), the
// partial leaves at the two ends of a range expanded by a document-array backend
// (plain DA, SA-Phi, RLCSA) instead of the grammar.
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

#include "dret/backend/backend.h"
#include "dret/backend/sa_phi.h"
#include "dret/collection.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/doc_list/sampled_tree_core.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/slp/slp_components.h"

namespace dret::gcda {

// The sampler of GCDA over a backend (sampled_tree_core.h): the cell's sampled tree stored alone
// (da_sampled_tree) and the document lists of its sampled nodes; the cover reads the tree, not the
// backend.
template <typename TStorage, uint8_t t_width, typename TSampledTree, typename TSets>
class SampledTree : public FieldsComponent<SampledTree<TStorage, t_width, TSampledTree, TSets>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<SampledTree, TStorage, t_width>;
  using typename Base::size_type;
  using Sets = TSets;

  static constexpr bool kBackendSharesCell = false;
  static constexpr const char* kBackendName = "get_docs";
  static constexpr const char* kBackendPrefix = "get_docs_";

  explicit SampledTree(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  SampledTree() = default;

  template <typename TBackend>
  void computeCoverFull(const TBackend& /*t_backend*/,
                        std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const {
    auto report = [&t_nodes](const auto& _value) { t_nodes.emplace_back(_value); };
    const auto range = grammar::ComputeCoverFromBottom(*tree_, t_sp, t_ep, report);
    CoverToRanges(t_sp, t_ep, range, !t_nodes.empty(), t_raw_ranges);
  }

  template <typename TReport>
  void forEachInNode(std::size_t t_i, TReport& t_report) const {
    for (const auto& d : (*sets_)[t_i]) t_report(static_cast<std::size_t>(d));
  }

  uint32_t block_size() const { return block_size_; }

  float storing_factor() const { return storing_factor_; }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto cell = [&t_self] { return SampledTreeCell{t_self.block_size_, t_self.storing_factor_}; };
    return std::tuple{
        Item("sampled_tree", t_self.tree_, [cell](const JSON& k) { return CellKey(k, conf::kDaSampledTree, cell()); }),
        Item("sets", t_self.sets_, [cell](const JSON& k) { return NodeDocListsKey<TSets>(k, cell()); })};
  }

  const TSampledTree* tree_ = nullptr;
  const TSets* sets_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

// The sampled tree and the node lists come from the cell's GCDA build (which reads the DA grammar).
template <typename TStorage, uint8_t t_width, typename TSampledTree, typename TSets>
void construct(SampledTree<TStorage, t_width, TSampledTree, TSets>& t_sampler, Config& t_config) {
  using namespace conf;
  const SampledTreeCell cell{t_sampler.block_size(), t_sampler.storing_factor()};
  if (!sdsl::cache_file_exists<TSampledTree>(CellKey(t_config.keys, kDaSampledTree, cell), t_config)) {
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
    LoadOrBuildSampledTree(t_config, filepath_da, cell.block_size, cell.storing_factor);
  }
  if (const auto key = NodeDocListsKey<TSets>(t_config.keys, cell); !sdsl::cache_file_exists<TSets>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    TSets sets;
    construct(sets, t_config, cell.block_size, cell.storing_factor);
  }
}

// TMergeSets is unused (the listing bitmap deduplicates); it stays for the existing spellings.
template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TBackend = backend::SAPhi<TStorage, TAlphabet::int_width>,
          typename TSets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                             true,
                                             grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TSampledTree = grammar::SampledSLP<>,
          typename TMergeSets = MergeSetsBinaryTreeFunctor>
using DocListIdxGCDABackend =
    DocListIdx<TStorage,
               TAlphabet,
               TCountIdx,
               SampledTreeCore<TStorage,
                               TAlphabet::int_width,
                               TBackend,
                               SampledTree<TStorage, TAlphabet::int_width, TSampledTree, TSets>>>;

}  // namespace dret::gcda
