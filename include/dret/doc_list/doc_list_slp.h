//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// Brute-GCDA: non-sampled GCDA index. Drives an SA-range search through a
// bare `grammar::SLP<>` (no sampled tree, no precomputed covers, no
// `GCChunks` doc-set storage), expanding the matched range and de-duping
// the resulting doc-id stream. It also takes a sampled GCDA grammar (cached,
// on-demand, compact) without its document lists: brute force over the
// grammar the RMQ/PDL GCDA backend loads.
//
// Compared to `dret::gcda::DocListIdxGCDA`, this class:
//   - Has no `block_size` / `storing_factor` knobs (those control the
//     sampled tree, which is absent here).
//   - Uses `grammar::ComputeSpanCover` to split `[sp, ep)` into a small
//     set of cover variables (each spans an exact subrange), then
//     expands each via `grammar::ExpandSLPForward` to enumerate
//     terminals = doc-ids.
//   - Reports each document the first time the expansion yields it, with
//     the bitmap every brute-force index shares (distinct_docs.h).
//
// `Search` cost is O(|cover| * log n) for the cover plus O((ep - sp) *
// height) for the expansions plus O(n_doc / 64) for the bitmap. There is no
// `getDocSet` shortcut because there are no precomputed sets — the whole
// [sp, ep) is walked.
//

#pragma once

#include <algorithm>
#include <filesystem>
#include <format>
#include <functional>
#include <iterator>
#include <vector>

#include <sdsl/construct_sa.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/memory_management.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "dret/collection.h"
#include "dret/construct_base.h"
#include "dret/doc_list/distinct_docs.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/index_base.h"
#include "dret/slp/bare_grammar.h"
#include "dret/slp/differential_slp.h"
#include "dret/slp/slp_components.h"

namespace dret {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>
class DocListIdxSLP : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using StorageBase = DocListIndexExtStorage<TStorage, TAlphabet>;
  using StorageBaseImpl = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using TPattern = typename TAlphabet::string_type;
  using TDocId = std::size_t;
  using size_type = std::size_t;

  // t_spacing is the sample spacing of a differential TSLP, or the block size of a
  // sampled one (with t_storing_factor); a plain one ignores both.
  explicit DocListIdxSLP(const TStorage& t_storage, uint32_t t_spacing = kDiffBlockSize,
                         float t_storing_factor = kNolistsStoringFactor)
      : StorageBase(t_storage), count_idx_(t_storage), spacing_(t_spacing), storing_factor_(t_storing_factor) {}

  DocListIdxSLP(const TStorage& t_storage, const TCountIdx& t_count_idx)
      : StorageBase(t_storage), count_idx_(t_count_idx) {}

  DocListIdxSLP() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = count_idx_.Count(t_pattern);
    if (sp >= ep) return;

    DistinctDocs seen(n_doc_);
    auto add = [&seen, &t_report](auto v) {
      if (seen.Insert(static_cast<std::size_t>(v))) t_report(static_cast<TDocId>(v));
    };

    // Plain SLP → span-cover expansion; differential SLP → range decompression;
    // sampled SLP → expansion through its sampled tree.
    SlpNsExpandRange(*slp_, sp, ep, add);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    sdsl::int_vector<64> n_doc_vec(1, n_doc_);
    return count_idx_.serialize(out, child, "count_idx") + slp_->serialize(out, child, "slp")
           + sdsl::serialize(n_doc_vec, out, child, "n_doc");
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    if (slp_) append(r, "slp", sdsl::size_in_bytes(*slp_));
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    return r;
  }

  const TCountIdx& count_idx() const { return count_idx_; }

  uint32_t spacing() const { return spacing_; }

  float storing_factor() const { return storing_factor_; }

 protected:
  void loadInner(typename StorageBaseImpl::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit([this](auto&& tt_source) { count_idx_.load(tt_source.get()); }, t_source);

    slp_ = this->template loadComponentsPtr<TSLP>(SlpNsTraits<TSLP>::components(t_keys, spacing_, storing_factor_),
                                                  t_source);

    n_doc_ = this->loadNDoc(t_keys, t_source);
  }

  TCountIdx count_idx_;
  const TSLP* slp_ = nullptr;
  uint32_t spacing_ = kDiffBlockSize;
  float storing_factor_ = kNolistsStoringFactor;
  std::size_t n_doc_ = 0;
};

template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TSLP>
void construct(DocListIdxSLP<TStorage, TAlphabet, TCountIdx, TSLP>& t_index, Config& t_config) {
  using namespace conf;

  EnsureDA<TAlphabet::int_width>(t_config);

  // Build the grammar's components — the plain CNF grammar, the base
  // differential SLP (bare-diff), or a sampled GCDA grammar — dispatched by TSLP.
  SlpNsTraits<TSLP>::build(t_config, t_index.spacing(), t_index.storing_factor());
  EnsureNDoc<TAlphabet::int_width>(t_config);

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  t_index.load(t_config);
}

}  // namespace dret
