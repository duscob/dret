//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/21/26.
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

#include "dret/collection.h"
#include "dret/slp/combined_slp_with_unit_cover.h"
#include "dret/slp/compact_bp_slp.h"
#include "dret/slp/compact_louds_slp.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_sampled_tree_base.h"
#include "dret/index_base.h"
#include "dret/slp/slp_components.h"
#include "dret/slp/gcda_grammar.h"
#include "dret/slp/slp_tools.h"

namespace dret {

namespace gcda {

class MergeSetsBinaryTreeFunctor;

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
class DocListIdxGCDA
    : public DLSampledTreeScheme<DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TSLP, TSLPSets, TMergeSets>,
                                 TMergeSets,
                                 typename TAlphabet::string_type>,
      public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  using SchemeBase = DLSampledTreeScheme<DocListIdxGCDA, TMergeSets, typename TAlphabet::string_type>;
  using StorageBase = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using typename SchemeBase::size_type;

  explicit DocListIdxGCDA(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_storage),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor) {}

  DocListIdxGCDA(const TStorage& t_storage,
                 const TCountIdx& t_count_idx,
                 uint32_t t_block_size = 512,
                 float t_storing_factor = 4)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_count_idx),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor) {}

  DocListIdxGCDA() = default;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return count_idx_.serialize(out, child, "count_idx") + slp_->serialize(out, child, "slp")
           + slp_sets_->serialize(out, child, "slp_docs");
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    if (slp_) collectSizes(r, *slp_, "slp_");
    if (slp_sets_) append(r, "slp_sets", sdsl::size_in_bytes(*slp_sets_));
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    return r;
  }

  const TCountIdx& count_idx() const {
    return count_idx_;
  }

  const uint32_t& block_size() const {
    return block_size_;
  }

  const float& storing_factor() const {
    return storing_factor_;
  }

 protected:
  // The steps of DLSampledTreeScheme::Query.
  friend SchemeBase;

  std::pair<std::size_t, std::size_t> count(const typename SchemeBase::TPattern& t_pattern) const {
    return count_idx_.Count(t_pattern);
  }

  void computeCoverFull(std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const {
    auto report = [&t_nodes](const auto& _value) { t_nodes.emplace_back(_value); };
    const auto range = grammar::ComputeCoverFromBottom(*slp_, t_sp, t_ep, report);
    SchemeBase::CoverToRanges(t_sp, t_ep, range, !t_nodes.empty(), t_raw_ranges);
  }

  // A reference to the sink, not the sink: some expansions take their report as const&.
  template <typename TSink>
  void getDocs(std::size_t t_sp, std::size_t t_ep, TSink& t_sink) const {
    auto position = [&t_sink](auto d) { t_sink.Position(static_cast<std::size_t>(d)); };
    ExpandSLP(*slp_, t_sp, t_ep, position);
  }

  // The node's list as the sets codec gives it (a grammar-compressed list expands into a vector; a
  // plain one is read in place).
  template <typename TSink>
  void forEachInNode(std::size_t t_i, TSink& t_sink) const {
    for (const auto& d : (*slp_sets_)[t_i]) t_sink.Listed(static_cast<std::size_t>(d));
  }

  void loadInner(typename StorageBase::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit(
        [this](auto&& tt_source) {
          count_idx_.load(tt_source.get());
        },
        t_source);

    const SampledTreeCell cell{block_size_, storing_factor_};
    slp_ = this->template loadComponentsPtr<TSLP>(CacheComponents(TSLP{}, t_keys, cell), t_source);
    slp_sets_ = this->template loadItemPtr<TSLPSets>(NodeDocListsKey<TSLPSets>(t_keys, cell), t_source, true);
  }

  TCountIdx count_idx_;
  const TSLP* slp_ = nullptr;
  const TSLPSets* slp_sets_ = nullptr;

  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TSLP,
          typename TSLPSets,
          typename TMergeSets>
void construct(DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TSLP, TSLPSets, TMergeSets>& t_index, Config& t_config) {
  using namespace conf;
  using Traits = GCDAVariantTraits<TSLP>;

  EnsureDA<TAlphabet::int_width>(t_config);

  const SampledTreeCell cell{t_index.block_size(), t_index.storing_factor()};

  if (const auto components = CacheComponents(TSLP{}, t_config.keys, cell); !ComponentsExist(components, t_config)) {
    auto event = sdsl::memory_monitor::event(components.back().key);
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
    TSLP slp;
    Traits::constructSLP(slp, t_config, filepath_da, t_index.block_size(), t_index.storing_factor());
  }

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  // The lists are the same for every representation of the cell (built while
  // sampling the tree, which every representation shares).
  if (const auto key_docs = NodeDocListsKey<TSLPSets>(t_config.keys, cell);
      !sdsl::cache_file_exists<TSLPSets>(key_docs, t_config)) {
    auto event = sdsl::memory_monitor::event(key_docs);
    TSLPSets slp_sets;
    construct(slp_sets, t_config, t_index.block_size(), t_index.storing_factor());
  }

  t_index.load(t_config);
}

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

// Backward-compatible alias. dret::dgcda::DocListIdxDGCDA locks TSLP to
// DifferentialLightSLP<> and forwards every other template parameter and
// default to gcda::DocListIdxGCDA. Existing consumers that wrote
// dret::dgcda::DocListIdxDGCDA<...> keep working unchanged; ADL on a
// DocListIdxGCDA<..., DifferentialLightSLP<>, ...> argument resolves to the
// unified gcda::construct(DocListIdxGCDA&, Config&) above.
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
