//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/16/26.
//
// DocListIdxPDL: PDL document listing, the sampled-tree core (sampled_tree_core.h)
// with PDLTree: a PDLTreeCore<..., TStoredSetCodec> for the sparse suffix-tree
// navigation and the per-slot stored sets, and a document-array backend for the
// ranges no stored node covers (default: backend::PlainDA; see backend/backend.h).
// The range comes from the r-index count (default: sri::RIndexCount), or from a
// backend that counts (the RLCSA, as in the original PDL; doc_list_index.h).
//
// The three PlainCodec / RPCodec / BCCodec variants are exposed as alias
// templates DocListIdxPDLPlain / DocListIdxPDLRP / DocListIdxPDLBC. Each codec
// names its lists component through a static constexpr kListsKey; the core is
// cached as components (see CacheComponents in pdl/tree_core.h).
//

#pragma once

#include <type_traits>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/construct_isa.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "sr-index/r_index.h"

#include "dret/backend/backend.h"
#include "dret/backend/da.h"
#include "dret/collection.h"
#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/doc_list/sampled_tree_core.h"
#include "dret/index_base.h"
#include "dret/pdl/build_pdl_core.h"
#include "dret/pdl/set_codecs.h"
#include "dret/pdl/storage_policy.h"
#include "dret/pdl/tree_builder.h"
#include "dret/pdl/tree_core.h"
#include "dret/size_report.h"

namespace dret::pdl {

// The sampler of PDL (sampled_tree_core.h): the sparse suffix tree of the collection, collapsed by the
// block size, with the stored sets of its selected nodes under a storage policy (PDLTreeCore, cached as
// components). Its cover can leave several uncovered ranges, between the selected nodes.
template <typename TStorage, uint8_t t_width, typename TCore>
class PDLTree : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;
  using Core = TCore;

  static constexpr bool kBackendSharesCell = false;
  static constexpr const char* kBackendName = "get_docs";
  static constexpr const char* kBackendPrefix = "get_docs_";

  explicit PDLTree(const TStorage& t_storage,
                   uint32_t t_block_size = 512,
                   float t_storing_factor = 4.0f,
                   StoragePolicy t_policy = StoragePolicy::OccurrenceWeighted)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor), policy_(t_policy) {}

  PDLTree() = default;

  // The selected nodes inside [sp, ep) and the gaps between them (PDLTreeCore).
  template <typename TBackend>
  void computeCoverFull(const TBackend& /*t_backend*/,
                        std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const {
    if (core_) core_->computeCoverFull(t_sp, t_ep, t_raw_ranges, t_nodes);
  }

  // The stored set, through the codec.
  template <typename TReport>
  void forEachInNode(std::size_t t_codec_slot, TReport& t_report) const {
    if (core_) core_->forEachInSet(t_codec_slot, t_report);
  }

  uint32_t block_size() const { return block_size_; }

  float storing_factor() const { return storing_factor_; }

  StoragePolicy policy() const { return policy_; }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return core_ ? core_->serialize(out, child, "core") : 0;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (core_)
      for (const auto& f : core_->GetSizeReport()) append(r, "core_" + f.name, f.bytes);
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    core_ = this->template loadComponentsPtr<TCore>(
        CacheComponents(TCore{}, t_keys, block_size_, storing_factor_, policy_), t_source);
  }

  const TCore* core_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4.0f;
  StoragePolicy policy_ = StoragePolicy::OccurrenceWeighted;
};

// The LCP array the sparse suffix tree is built from (the DA comes from the core's construct), then
// the tree core.
template <typename TStorage, uint8_t t_width, typename TCore>
void construct(PDLTree<TStorage, t_width, TCore>& t_sampler, Config& t_config) {
  using namespace conf;

  EnsureLCP<t_width>(t_config);

  // PDLTreeCore as components: the tree under the block size, the selected
  // nodes and their lists under the block size, the policy and, for a policy
  // that reads it, the storing factor.
  const auto components =
      CacheComponents(TCore{}, t_config.keys, t_sampler.block_size(), t_sampler.storing_factor(), t_sampler.policy());

  if (!ComponentsExist(components, t_config)) {
    auto event = sdsl::memory_monitor::event(components.back().key);
    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);
    sdsl::int_vector<> lcp;
    sdsl::load_from_cache(lcp, sdsl::conf::KEY_LCP, t_config);

    const std::size_t n = da.size();
    const std::size_t n_doc = ReadNDoc<t_width>(t_config);

    // PDL construction pipeline (pdl/tree_builder.h) plus the bridge into the
    // compact PDLTreeCore (build_pdl_core.h).
    BuilderPool pool;
    // BuildSparseSuffixTree iterates i in [1, n] and reads lcp(n) as a
    // flush sentinel; sdsl::construct_lcp_kasai writes only n entries, so
    // wrap with a bounds-safe accessor returning 0 past the end.
    auto lcp_fn = [&lcp](std::size_t i) -> std::size_t {
      return i < lcp.size() ? static_cast<std::size_t>(lcp[i]) : 0;
    };
    auto* root = BuildSparseSuffixTree(pool, lcp_fn, n);
    CollapseSubtreesByBlockSize(root, t_sampler.block_size());
    // Reclaim the interval nodes just detached by the collapse before the
    // memory-heavy doc-set phase (the arena keeps all nodes alive otherwise).
    pool.retain(root);
    InsertExplicitLeaves(pool, root);
    ComputeDocSetsBottomUp(
        root, n_doc, [&da](std::size_t i) { return static_cast<std::size_t>(da[i]); });
    ApplyStoragePolicy(root, t_sampler.policy(), t_sampler.storing_factor(), n_doc);
    const std::size_t n_nodes = AssignNodeIds(root);

    TCore core;
    BuildPDLTreeCoreFromBuilder(
        core, root, n_nodes, n_doc, t_sampler.block_size(), t_sampler.storing_factor(), t_sampler.policy());
    StoreComponents(
        core,
        CacheComponents(core, t_config.keys, t_sampler.block_size(), t_sampler.storing_factor(), t_sampler.policy()),
        t_config);
  }
}

// PDL: the sampled-tree core with PDLTree, parameterised on the stored-set codec. TMergeSets is unused
// (the listing bitmap deduplicates); it stays for the existing spellings.
template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TBackend = backend::PlainDA<TStorage, TAlphabet::int_width>,
          typename TStoredSetCodec = PlainCodec<>,
          typename TBitvector = sdsl::sd_vector<>,
          typename TIntVector = sdsl::int_vector<>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
using DocListIdxPDL = DocListIdx<TStorage,
                                 TAlphabet,
                                 TCountIdx,
                                 SampledTreeCore<TStorage,
                                                 TAlphabet::int_width,
                                                 TBackend,
                                                 PDLTree<TStorage,
                                                         TAlphabet::int_width,
                                                         PDLTreeCore<TBitvector,
                                                                     typename TBitvector::rank_1_type,
                                                                     typename TBitvector::select_1_type,
                                                                     TIntVector,
                                                                     TStoredSetCodec>>>>;

// The three codecs: plain, Re-Pair and BC lists.

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TBackend = backend::PlainDA<TStorage, TAlphabet::int_width>,
          typename TBitvector = sdsl::sd_vector<>,
          typename TIntVector = sdsl::int_vector<>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
using DocListIdxPDLPlain =
    DocListIdxPDL<TStorage, TAlphabet, TCountIdx, TBackend, PlainCodec<>, TBitvector, TIntVector, TMergeSets>;

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TBackend = backend::PlainDA<TStorage, TAlphabet::int_width>,
          typename TBitvector = sdsl::sd_vector<>,
          typename TIntVector = sdsl::int_vector<>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
using DocListIdxPDLRP =
    DocListIdxPDL<TStorage, TAlphabet, TCountIdx, TBackend, RPCodec<>, TBitvector, TIntVector, TMergeSets>;

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TBackend = backend::PlainDA<TStorage, TAlphabet::int_width>,
          typename TBitvector = sdsl::sd_vector<>,
          typename TIntVector = sdsl::int_vector<>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
using DocListIdxPDLBC =
    DocListIdxPDL<TStorage, TAlphabet, TCountIdx, TBackend, BCCodec<>, TBitvector, TIntVector, TMergeSets>;

}  // namespace dret::pdl
