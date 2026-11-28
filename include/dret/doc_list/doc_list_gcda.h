//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/21/26.
//
// GCDA: the sampled-tree core (sampled_tree_core.h) over the grammar-compressed DA. The grammar is the
// backend that expands the uncovered ends of a range (backend::GCDAGrammar), and its sampled parse tree
// gives the cover (GrammarTree); the stored lists of the sampled nodes are a separate component, shared
// by every representation of the cell. GCDA over another backend keeps the tree alone
// (doc_list_gcda_backend.h).
//

#pragma once

#include <stdexcept>
#include <type_traits>

#include <sdsl/construct_sa.hpp>

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "sr-index/r_index.h"

#include "dret/backend/grammar.h"
#include "dret/collection.h"
#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/doc_list/sampled_tree_core.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/slp/combined_slp_with_unit_cover.h"
#include "dret/slp/compact_bp_slp.h"
#include "dret/slp/compact_louds_slp.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/slp/gcda_grammar.h"
#include "dret/slp/slp_components.h"
#include "dret/slp/slp_tools.h"

namespace dret {

namespace gcda {

class MergeSetsBinaryTreeFunctor;

// The sampler of GCDA (sampled_tree_core.h): the sampled parse tree of the DA grammar, read from the
// backend's grammar (backend::GCDAGrammar, which also expands the uncovered ends), and the document
// lists of its sampled nodes. The grammar is the backend, so it is built for the sampler's cell.
template <typename TStorage, uint8_t t_width, typename TSets>
class GrammarTree : public FieldsComponent<GrammarTree<TStorage, t_width, TSets>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<GrammarTree, TStorage, t_width>;
  using typename Base::size_type;
  using Sets = TSets;

  static constexpr bool kBackendSharesCell = true;
  static constexpr const char* kBackendName = "slp";
  static constexpr const char* kBackendPrefix = "";

  explicit GrammarTree(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  GrammarTree() = default;

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

  // The node's list as the sets codec gives it (a grammar-compressed list expands into a vector; a
  // plain one is read in place).
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
    auto key = [&t_self](const JSON& k) {
      return NodeDocListsKey<TSets>(k, {t_self.block_size_, t_self.storing_factor_});
    };
    return std::tuple{Item("slp_sets", t_self.sets_, key)};
  }

  const TSets* sets_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

// The lists are the same for every representation of the cell (built while sampling the tree, which
// every representation shares).
template <typename TStorage, uint8_t t_width, typename TSets>
void construct(GrammarTree<TStorage, t_width, TSets>& t_sampler, Config& t_config) {
  const SampledTreeCell cell{t_sampler.block_size(), t_sampler.storing_factor()};
  if (const auto key_docs = NodeDocListsKey<TSets>(t_config.keys, cell);
      !sdsl::cache_file_exists<TSets>(key_docs, t_config)) {
    auto event = sdsl::memory_monitor::event(key_docs);
    TSets sets;
    construct(sets, t_config, cell.block_size, cell.storing_factor);
  }
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
  // The lists carry counts (sampled_tree_core.h).
  static constexpr bool kWeighted = true;

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

  // t_report(d, f) for every document d of the node's list, with its frequency f in the node.
  template <typename TReport>
  void forEachWeightedInNode(std::size_t t_i, TReport& t_report) const {
    const auto docs = (*sets_)[t_i];
    auto f = (*counts_)[t_i].begin();
    for (auto d = docs.begin(); d != docs.end(); ++d, ++f)
      t_report(static_cast<std::size_t>(*d), static_cast<std::size_t>(*f));
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
  GrammarTree<TStorage, t_width, TSets> lists(t_sampler.storage(), t_sampler.block_size(), t_sampler.storing_factor());
  construct(lists, t_config);
  const SampledTreeCell cell{t_sampler.block_size(), t_sampler.storing_factor()};
  if (!sdsl::cache_file_exists<NodeDocCounts>(NodeDocCountsKey(t_config.keys, cell), t_config)) {
    const auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[conf::kDA].get<std::string>(), t_config);
    NodeDocCounts counts;
    construct(counts, t_config, filepath_da, cell.block_size, cell.storing_factor);
  }
}

// GCDA: the sampled-tree core over the grammar, with GrammarTree. TMergeSets is unused (the listing
// bitmap deduplicates, distinct_docs.h); it stays for the existing spellings of the type.
template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TSLPSets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                                true,
                                                grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TMergeSets = MergeSetsBinaryTreeFunctor>
using DocListIdxGCDA = DocListIdx<TStorage,
                                  TAlphabet,
                                  TCountIdx,
                                  SampledTreeCore<TStorage,
                                                  TAlphabet::int_width,
                                                  backend::GCDAGrammar<TStorage, TAlphabet::int_width, TSLP>,
                                                  GrammarTree<TStorage, TAlphabet::int_width, TSLPSets>>>;

// GCDA with the counts of its lists: computes document frequencies (Full-GCDA).
template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TSLPSets = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>
using DocFreqIdxGCDA = DocListIdx<TStorage,
                                  TAlphabet,
                                  TCountIdx,
                                  SampledTreeCore<TStorage,
                                                  TAlphabet::int_width,
                                                  backend::GCDAGrammar<TStorage, TAlphabet::int_width, TSLP>,
                                                  CountedGrammarTree<TStorage, TAlphabet::int_width, TSLPSets>>>;

//~~~~~~~


class MergeSetsBinaryTreeFunctor {
 public:
  MergeSetsBinaryTreeFunctor() = default;

  template <typename TII, typename TSets, typename TResult>
  void operator()(TII _first, TII _last, const TSets& _sets, TResult& _result) const {
    auto default_set_union = [](auto _first1, auto _last1, auto _first2, auto _last2, auto _result) -> auto {
      return std::set_union(_first1, _last1, _first2, _last2, _result);
    };

    (*this)(_first, _last, _sets, _result, default_set_union);
  }

  template <typename _II, typename _Sets, typename _Result, typename _SetUnion>
  void operator()(_II _first, _II _last, const _Sets& _sets, _Result& _result, const _SetUnion& _set_union) const {
    _Result tmp_merge;

    auto merge_tmp = [&tmp_merge, &_set_union](const auto& set1, const auto& set2) {
      tmp_merge.resize(set1.size() + set2.size());
      auto last_it = _set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), tmp_merge.begin());
      tmp_merge.resize(last_it - tmp_merge.begin());
    };

    auto length = std::distance(_first, _last);
    if (length == 1) {
      merge_tmp(_result, _sets(*_first));

      _result.swap(tmp_merge);

      return;
    }

    std::vector<std::pair<uint8_t, _Result>> part_results = {{1, {}}};
    part_results.front().second.swap(_result);

    while (part_results.size() != 1 || _first != _last) {
      std::size_t size;
      while ((size = part_results.size()) > 1
             && (part_results[size - 1].first == part_results[size - 2].first || _first == _last)) {
        merge_tmp(part_results[size - 1].second, part_results[size - 2].second);

        part_results[size - 2].second.swap(tmp_merge);
        ++part_results[size - 2].first;
        part_results.pop_back();
      }

      if (_first != _last) {
        auto next = _first + 1;
        if (next == _last) {
          merge_tmp(part_results.back().second, _sets(*_first));

          part_results.back().second.swap(tmp_merge);
          ++_first;
        } else {
          merge_tmp(_sets(*_first), _sets(*next));

          part_results.emplace_back(1, std::move(tmp_merge));
          _first += 2;
        }
      }
    }

    _result.swap(part_results.front().second);
  }
};

}  // namespace gcda

//~~~~~~~

// GCDA over the differential grammar: DocListIdxGCDA with TSLP = DifferentialLightSLP<>.
namespace dgcda {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = DifferentialLightSLP<>,
          typename TSLPSets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                                  true,
                                                  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
using DocListIdxDGCDA = gcda::DocListIdxGCDA<TStorage, TAlphabet, TCountIdx,
                                              TSLP, TSLPSets, TMergeSets>;

}  // namespace dgcda

}  // namespace dret
