//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// DocListIdxPDLBC: PDL document-listing index whose stored sets use the
// BCCodec (biClique-coded blocks built via Cecilia Hernandez's
// vnmextract dense-subgraph extractor, see external/dsextract/). Mirrors
// DocListIdxPDLPlain / DocListIdxPDLRP in shape — only the stored-set
// codec, the JSON key path (kBC instead of kPlain/kRP), and the
// cache-key prefix differ. Search() is inherited from
// DLSampledTreeScheme.
//
// Task 25 of docs/pdl_indexes_tasks.md.
//

#pragma once

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

#include "../config.h"
#include "../construct_base.h"
#include "../doc_list_sampled_tree.h"
#include "../doc_list_sampled_tree_gcda.h"
#include "../index_base.h"
#include "../size_report.h"
#include "build_pdl_core.h"
#include "doc_list_pdl_plain.h"  // for ReadNDocFromText
#include "get_docs.h"
#include "set_codecs.h"
#include "storage_policy.h"
#include "tree_builder.h"
#include "tree_core.h"

namespace dret::pdl {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TGetDocs = PDLGetDocsDA<TStorage, TAlphabet::int_width>,
          typename TStoredSetCodec = BCCodec<>,
          typename TBitvector = sdsl::sd_vector<>,
          typename TIntVector = sdsl::int_vector<>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
class DocListIdxPDLBC
    : public DLSampledTreeScheme<TMergeSets, typename TAlphabet::string_type>,
      public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  using SchemeBase = DLSampledTreeScheme<TMergeSets, typename TAlphabet::string_type>;
  using StorageBase = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using TCore = PDLTreeCore<TBitvector,
                            typename TBitvector::rank_1_type,
                            typename TBitvector::select_1_type,
                            TIntVector,
                            TStoredSetCodec>;
  using typename SchemeBase::size_type;

  DocListIdxPDLBC() = default;

  explicit DocListIdxPDLBC(const TStorage& t_storage,
                           uint32_t t_block_size = 512,
                           float t_storing_factor = 4.0f,
                           StoragePolicy t_policy = StoragePolicy::OccurrenceWeighted)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_storage),
        get_docs_(typename TGetDocs::Inner{t_storage}),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor),
        policy_(t_policy) {}

  DocListIdxPDLBC(const TStorage& t_storage,
                  const TCountIdx& t_count_idx,
                  uint32_t t_block_size = 512,
                  float t_storing_factor = 4.0f,
                  StoragePolicy t_policy = StoragePolicy::OccurrenceWeighted)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_count_idx),
        get_docs_(typename TGetDocs::Inner{t_storage}),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor),
        policy_(t_policy) {}

  size_type serialize(std::ostream& out,
                      sdsl::structure_tree_node* v,
                      const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    std::size_t bytes = 0;
    bytes += count_idx_.serialize(out, child, "count_idx");
    bytes += get_docs_.inner().serialize(out, child, "get_docs");
    if (core_) bytes += core_->serialize(out, child, "core");
    return bytes;
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    auto gd = get_docs_.inner().GetSizeReport();
    for (const auto& f : gd) append(r, "get_docs_" + f.name, f.bytes);
    if (core_) {
      auto core_r = core_->GetSizeReport();
      for (const auto& f : core_r) append(r, "core_" + f.name, f.bytes);
    }
    return r;
  }

  const TCountIdx& count_idx() const { return count_idx_; }
  uint32_t block_size() const { return block_size_; }
  float storing_factor() const { return storing_factor_; }
  StoragePolicy policy() const { return policy_; }

 protected:
  std::pair<std::size_t, std::size_t> count(const typename SchemeBase::TPattern& t_pattern) const override {
    return count_idx_.Count(t_pattern);
  }

  std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> computeCover(
      std::size_t t_sp,
      std::size_t /*t_ep*/) const override {
    return {{t_sp, t_sp}, {}};
  }

  void computeCoverFull(std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const override {
    if (core_) core_->computeCoverFull(t_sp, t_ep, t_raw_ranges, t_nodes);
  }

  void getDocs(std::size_t t_sp,
               std::size_t t_ep,
               const std::function<void(typename SchemeBase::TDocId)>& t_report) const override {
    auto wrapper = [&t_report](std::size_t d) { t_report(d); };
    get_docs_.getDocs(t_sp, t_ep, wrapper);
  }

  std::vector<typename SchemeBase::TDocId> getDocSet(std::size_t t_codec_slot) const override {
    if (!core_) return {};
    return core_->getDocSet(t_codec_slot);
  }

  void loadInner(typename StorageBase::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit(
        [this](auto&& tt_source) {
          count_idx_.load(tt_source.get());
          get_docs_.inner().load(tt_source.get());
        },
        t_source);

    auto key_prefix = std::format("{}-{}_pdl_bc_{}_", block_size_, storing_factor_,
                                  static_cast<int>(policy_));
    auto key_core = key_prefix + t_keys[kPDL][kBC][kSets].get<std::string>();
    core_ = this->template loadItemPtr<TCore>(key_core, t_source, true);
  }

  TCountIdx count_idx_;
  TGetDocs get_docs_;
  const TCore* core_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4.0f;
  StoragePolicy policy_ = StoragePolicy::OccurrenceWeighted;
};

template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TGetDocs,
          typename TStoredSetCodec,
          typename TBitvector,
          typename TIntVector,
          typename TMergeSets>
void construct(DocListIdxPDLBC<TStorage, TAlphabet, TCountIdx, TGetDocs,
                               TStoredSetCodec, TBitvector, TIntVector, TMergeSets>& t_index,
               Config& t_config) {
  using namespace conf;
  using Index = std::remove_reference_t<decltype(t_index)>;
  using TCore = typename Index::TCore;
  constexpr uint8_t t_width = TAlphabet::int_width;

  if (!cache_file_exists(t_config.keys[kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<t_width>(t_config);
  }
  if (const auto key = t_config.keys[kSA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    sdsl::construct_sa<t_width>(t_config);
  }
  if (const auto key = t_config.keys[kDocEnds].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocEnd<t_width, sdsl::sd_vector<>>(t_config);
  }
  if (const auto key = t_config.keys[kDA].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray(t_config);
  }
  if (!cache_file_exists(sdsl::conf::KEY_LCP, t_config)) {
    auto event = sdsl::memory_monitor::event(std::string(sdsl::conf::KEY_LCP));
    if (!cache_file_exists(sdsl::conf::KEY_ISA, t_config)) {
      sdsl::construct_isa(t_config);
    }
    sdsl::construct_lcp_kasai<t_width>(t_config);
  }

  const auto key_prefix = std::format("{}-{}_pdl_bc_{}_",
                                      t_index.block_size(),
                                      t_index.storing_factor(),
                                      static_cast<int>(t_index.policy()));
  const auto key_core = key_prefix + t_config.keys[kPDL][kBC][kSets].get<std::string>();

  if (!sdsl::cache_file_exists<TCore>(key_core, t_config)) {
    auto event = sdsl::memory_monitor::event(key_core);
    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);
    sdsl::int_vector<> lcp;
    sdsl::load_from_cache(lcp, sdsl::conf::KEY_LCP, t_config);

    const std::size_t n = da.size();
    const std::size_t n_doc = ReadNDocFromText(t_config);

    BuilderPool pool;
    auto lcp_fn = [&lcp](std::size_t i) -> std::size_t {
      return i < lcp.size() ? static_cast<std::size_t>(lcp[i]) : 0;
    };
    auto* root = BuildSparseSuffixTree(pool, lcp_fn, n);
    CollapseSubtreesByBlockSize(root, t_index.block_size());
    InsertExplicitLeaves(pool, root);
    ComputeDocSetsBottomUp(
        root, n_doc, [&da](std::size_t i) { return static_cast<std::size_t>(da[i]); });
    ApplyStoragePolicy(root, t_index.policy(), t_index.storing_factor());
    const std::size_t n_nodes = AssignNodeIds(root);

    TCore core;
    BuildPDLTreeCoreFromBuilder(core, root, n_nodes, n_doc,
                                t_index.block_size(), t_index.storing_factor(),
                                t_index.policy());
    sdsl::store_to_cache(core, key_core, t_config, true);
  }

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  t_index.load(t_config);
}

}  // namespace dret::pdl
