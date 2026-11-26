//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// GCDA with the counts of its node lists: how often each document of a sampled node's list occurs in
// the node, stored next to the lists (keyed by the cell), so the sampled-tree core can hand a
// frequency sink the nodes with their counts (CountedGrammarTree::ReportNode). The grammar, the tree
// and the lists are GCDA's (doc_list/doc_list_gcda.h, slp/gcda_grammar.h).
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>
#include <grammar/slp_metadata.h>

#include "dret/config.h"
#include "dret/doc_freq/frequencies.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/sampled_tree_core.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/slp/gcda_grammar.h"
#include "dret/slp/slp_components.h"

namespace dret::freq {

using gcda::LoadRawNodeDocLists;
using gcda::NodeDocListsKey;

// The frequencies of the documents of each sampled node's list: how often each occurs in the node's
// expansion, in the order of the list (GCDA's counts, for document frequencies). Bit-packed.
using NodeDocCounts = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;

inline std::string NodeDocCountsKey(const JSON& t_keys, const SampledTreeCell& t_cell) {
  return CellKey(t_keys, conf::kDaNodeDocCounts, t_cell);
}

namespace internal {

// grammar::AddSet (grammar/sampled_slp.h), which inserts the sorted distinct documents of a node's
// expansion, that also inserts how often each occurs.
class AddCountedSet {
 public:
  AddCountedSet(grammar::Chunks<>& t_docs, grammar::Chunks<>& t_counts) : docs_(t_docs), counts_(t_counts) {}

  template <typename TSLP>
  void operator()(const TSLP& t_slp, std::size_t t_var) {
    auto span = t_slp.Span(t_var);
    std::sort(span.begin(), span.end());
    std::vector<uint32_t> docs, counts;
    for (std::size_t i = 0; i < span.size();) {
      std::size_t j = i;
      while (j < span.size() && span[j] == span[i]) ++j;
      docs.push_back(span[i]);
      counts.push_back(j - i);
      i = j;
    }
    docs_.Insert(docs.begin(), docs.end());
    counts_.Insert(counts.begin(), counts.end());
  }

  template <typename TSLP, typename TNodes, typename TRanges>
  void operator()(const TSLP& t_slp, std::size_t t_var, const TNodes&, std::size_t, const TRanges&, const TRanges&) {
    (*this)(t_slp, t_var);
  }

 private:
  grammar::Chunks<>& docs_;
  grammar::Chunks<>& counts_;
};

}  // namespace internal

// The cell's sampled tree again, recording the counts as its node lists are computed. The tree is a
// function of the grammar and the cell, so it has the same nodes in the same order as the build that
// cached the lists; the lists are compared with the cached ones to make sure.
inline void construct(NodeDocCounts& t_counts,
                      Config& t_config,
                      const std::string& t_datafile,
                      uint32_t t_block_size,
                      float t_storing_factor) {
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto key = NodeDocCountsKey(t_config.keys, cell);
  auto event = sdsl::memory_monitor::event(key);

  grammar::CombinedSLP<> cslp(LoadOrBuildDaCnfGrammar(t_config, t_datafile));
  grammar::Chunks<> docs, counts;
  internal::AddCountedSet add(docs, counts);
  cslp.Compute(t_block_size,
               add,
               add,
               grammar::MustBeSampled<decltype(docs)>(grammar::AreChildrenTooBig(docs, t_storing_factor)));
  if (!(docs == LoadRawNodeDocLists(t_config, t_block_size, t_storing_factor)))
    throw std::runtime_error("GCDA document counts: the rebuilt tree's lists differ from the cached ones (" + key + ")");

  auto bit_compress = [](sdsl::int_vector<>& v) { sdsl::util::bit_compress(v); };
  t_counts = NodeDocCounts(counts, bit_compress, bit_compress);
  sdsl::store_to_cache(t_counts, key, t_config, true);
}

// GrammarTree with the frequencies of each list's documents (NodeDocCounts): a GCDA that computes
// document frequencies, Full-GCDA in the 2020 frequency code. The counts are a component of their own,
// keyed by the cell; the lists and the grammar are GCDA's.
template <typename TStorage, uint8_t t_width, typename TSets>
class CountedGrammarTree : public FieldsComponent<CountedGrammarTree<TStorage, t_width, TSets>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<CountedGrammarTree, TStorage, t_width>;
  using typename Base::size_type;
  using Sets = TSets;

  static constexpr bool kBackendSharesCell = true;
  static constexpr const char* kBackendName = "slp";
  static constexpr const char* kBackendPrefix = "";

  explicit CountedGrammarTree(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  CountedGrammarTree() = default;

  template <typename TBackend>
  void computeCoverFull(const TBackend& t_backend,
                        std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const {
    auto report = [&t_nodes](const auto& _value) { t_nodes.emplace_back(_value); };
    const auto range = grammar::ComputeCoverFromBottom(t_backend.slp(), t_sp, t_ep, report);
    CoverToRanges(t_sp, t_ep, range, !t_nodes.empty(), t_raw_ranges);
  }

  template <typename TReport>
  void forEachInNode(std::size_t t_i, TReport& t_report) const {
    for (const auto& d : (*sets_)[t_i]) t_report(static_cast<std::size_t>(d));
  }

  // The node's list into the sink (sampled_tree_core.h): every document with its frequency in the node
  // when the sink takes counts (Weighted), else every document.
  template <typename TSink>
  void ReportNode(std::size_t t_i, TSink& t_sink) const {
    if constexpr (requires { t_sink.Weighted(std::size_t{}, std::size_t{}); }) {
      const auto docs = (*sets_)[t_i];
      auto f = (*counts_)[t_i].begin();
      for (auto d = docs.begin(); d != docs.end(); ++d, ++f)
        t_sink.Weighted(static_cast<std::size_t>(*d), static_cast<std::size_t>(*f));
    } else {
      for (const auto& d : (*sets_)[t_i]) t_sink.Listed(static_cast<std::size_t>(d));
    }
  }

  uint32_t block_size() const { return block_size_; }

  float storing_factor() const { return storing_factor_; }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto cell = [&t_self] { return SampledTreeCell{t_self.block_size_, t_self.storing_factor_}; };
    return std::tuple{Item("slp_sets", t_self.sets_, [cell](const JSON& k) { return NodeDocListsKey<TSets>(k, cell()); }),
                      Item("slp_counts", t_self.counts_, [cell](const JSON& k) { return NodeDocCountsKey(k, cell()); })};
  }

  const TSets* sets_ = nullptr;
  const NodeDocCounts* counts_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

// GrammarTree's lists, then the counts, built from the grammar the lists were built from.
template <typename TStorage, uint8_t t_width, typename TSets>
void construct(CountedGrammarTree<TStorage, t_width, TSets>& t_sampler, Config& t_config) {
  gcda::GrammarTree<TStorage, t_width, TSets> lists(t_sampler.storage(), t_sampler.block_size(), t_sampler.storing_factor());
  construct(lists, t_config);
  const SampledTreeCell cell{t_sampler.block_size(), t_sampler.storing_factor()};
  if (!sdsl::cache_file_exists<NodeDocCounts>(NodeDocCountsKey(t_config.keys, cell), t_config)) {
    const auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[conf::kDA].get<std::string>(), t_config);
    NodeDocCounts counts;
    construct(counts, t_config, filepath_da, cell.block_size, cell.storing_factor);
  }
}

// GCDA with the counts of its lists: computes document frequencies (Full-GCDA in the 2020 code).
template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TSLPSets = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>
using DocFreqIdxGCDA = DocIndex<TStorage,
                                  TAlphabet,
                                  TCountIdx,
                                  SampledTreeCore<TStorage,
                                                  TAlphabet::int_width,
                                                  backend::GCDAGrammar<TStorage, TAlphabet::int_width, TSLP>,
                                                  CountedGrammarTree<TStorage, TAlphabet::int_width, TSLPSets>>>;

// The sampled-tree core with counted lists serves frequencies: the uncovered ends go position by
// position, the nodes with their counts.
template <typename TStorage, uint8_t t_width, typename TBackend, typename TStorage2, uint8_t t_width2, typename TSets>
struct ComputesFrequencies<SampledTreeCore<TStorage, t_width, TBackend, CountedGrammarTree<TStorage2, t_width2, TSets>>>
    : std::true_type {};

}  // namespace dret::freq
