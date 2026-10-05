//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// PDL sparse-suffix-tree core. Holds the compact navigation
// representation (per-node intervals + linked-list children encoded as
// id arrays + selected-node bitvector) and answers multi-range cover
// queries from DLSampledTreeScheme. Built either from a BuilderNode tree
// at construction (Task 28's construct() via build_pdl_core.h) or from
// disk via load() (Task 13).
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <stdexcept>
#include <string_view>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "dret/cache_components.h"
#include "dret/config.h"
#include "dret/size_report.h"
#include "dret/pdl/storage_policy.h"

namespace dret::pdl {

// Placeholder codec satisfying the minimum API used by PDLTreeCore. Replaced
// by Plain / RP / BC codecs (set_codecs.h, Tasks 15–17). Includes a no-op
// Expand so PDLTreeCore::getDocSet still instantiates against this default
// codec — querying then yields no documents per slot, which matches the
// "skeleton, no real data" semantics.
struct NullCodec {
  // Emits nothing, which is trivially sorted.
  static constexpr bool kExpandsSorted = true;

  template <typename TGetSetAt>
  void Build(std::size_t /*t_n_slots*/, TGetSetAt&& /*t_get_set_at*/,
             std::size_t /*t_n_doc*/) {}

  template <typename TReport>
  void Expand(std::size_t /*t_slot*/, std::size_t /*t_n_doc*/, TReport&& /*t_report*/) const {}

  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    sdsl::structure_tree::add_size(child, 0);
    return 0;
  }
  void load(std::istream& /*in*/) {}

  SizeReport GetSizeReport() const { return {}; }
};

template <typename TBitvector = sdsl::sd_vector<>,
          typename TBvRank = typename TBitvector::rank_1_type,
          typename TBvSelect = typename TBitvector::select_1_type,
          typename TIntVector = sdsl::int_vector<>,
          typename TStoredSetCodec = NullCodec>
class PDLTreeCore {
  template <typename, typename, typename, typename, typename>
  friend class PDLTreeCore;

 public:
  using TDocId = std::size_t;
  using size_type = std::size_t;

  PDLTreeCore() = default;

  // A core whose stored sets use another codec type -- in practice one cached
  // with RPCodecWithLengths: the tree is copied as is and only the codec is
  // converted, through TStoredSetCodec's converting constructor.
  template <typename TOtherCodec>
  explicit PDLTreeCore(const PDLTreeCore<TBitvector, TBvRank, TBvSelect, TIntVector, TOtherCodec>& t_other) {
    Assemble(t_other.node_starts_, t_other.node_ends_, t_other.first_child_, t_other.next_sibling_,
             t_other.selected_marker_, TStoredSetCodec(t_other.stored_sets_), t_other.n_doc_,
             t_other.block_size_, t_other.storing_factor_, t_other.policy_);
  }

  // Single-range cover; matches DLSampledTreeScheme::computeCover. Stub
  // remains for Tasks 23-25 that may need the legacy single-range API
  // alongside the multi-range computeCoverFull below.
  std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> computeCover(
      std::size_t t_sp,
      std::size_t /*t_ep*/) const {
    return {{t_sp, t_sp}, {}};
  }

  // Multi-range cover. Iterative DFS from the root: a node fully inside
  // [sp, ep) contributes either its codec slot (when selected) or — if
  // it's a leaf — its raw range; navigation-only internal nodes recurse
  // into their children. Partial overlaps always recurse (or yield the
  // overlap as raw if the node is a leaf).
  //
  // Output `t_nodes` carries codec slots (selected-rank of the node id),
  // not raw node ids — this matches DLSampledTreeScheme's getDocSet
  // contract, which indexes the stored-set codec.
  void computeCoverFull(std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const {
    if (n_nodes_ == 0 || t_sp >= t_ep) return;

    // Root id is n_nodes_ - 1 by AssignNodeIds' post-order convention.
    std::vector<std::size_t> stack{n_nodes_ - 1};
    while (!stack.empty()) {
      const std::size_t id = stack.back();
      stack.pop_back();
      const std::size_t n_sp = node_starts_[id];
      const std::size_t n_ep = node_ends_[id];
      if (n_ep <= t_sp || n_sp >= t_ep) continue;  // disjoint

      const std::size_t fc = first_child_[id];
      const bool is_leaf = (fc == kSentinel());
      const bool fully_inside = (t_sp <= n_sp && n_ep <= t_ep);

      if (fully_inside && selected_marker_[id]) {
        t_nodes.push_back(selected_rank_(id));
        continue;
      }
      if (is_leaf) {
        const std::size_t b = fully_inside ? n_sp : std::max(t_sp, n_sp);
        const std::size_t e = fully_inside ? n_ep : std::min(t_ep, n_ep);
        t_raw_ranges.emplace_back(b, e);
        continue;
      }
      // Navigation node OR partial-overlap internal: descend.
      for (std::size_t cid = fc; cid != kSentinel(); cid = next_sibling_[cid]) {
        stack.push_back(cid);
      }
    }
  }

  std::vector<TDocId> getDocSet(std::size_t t_codec_slot) const {
    std::vector<TDocId> out;
    out.reserve(n_doc_);
    stored_sets_.Expand(t_codec_slot, n_doc_,
                        [&out](std::size_t d) { out.push_back(d); });
    // DLSampledTreeScheme::Search feeds this straight into
    // MergeSetsBinaryTreeFunctor, which merges with std::set_union -- so the
    // result has to be sorted and unique. Codecs that already emit in
    // ascending order (Plain, RP) declare kExpandsSorted and pay nothing here;
    // BC concatenates per-rule runs and does not, so it is sorted once, at the
    // only place where the ordering contract is actually consumed.
    if constexpr (!TStoredSetCodec::kExpandsSorted) {
      std::sort(out.begin(), out.end());
      out.erase(std::unique(out.begin(), out.end()), out.end());
    }
    return out;
  }

  // Internal — populate the compact representation in one shot. Used by
  // build_pdl_core.h (construction) and load() (disk path). The rank /
  // select supports MUST be re-pointed at the new bitvector after every
  // assignment; doing it here keeps callers from forgetting.
  void Assemble(TIntVector t_node_starts,
                TIntVector t_node_ends,
                TIntVector t_first_child,
                TIntVector t_next_sibling,
                TBitvector t_selected_marker,
                TStoredSetCodec t_stored_sets,
                std::size_t t_n_doc,
                uint32_t t_block_size,
                float t_storing_factor,
                StoragePolicy t_policy) {
    node_starts_ = std::move(t_node_starts);
    node_ends_ = std::move(t_node_ends);
    first_child_ = std::move(t_first_child);
    next_sibling_ = std::move(t_next_sibling);
    selected_marker_ = std::move(t_selected_marker);
    stored_sets_ = std::move(t_stored_sets);
    n_doc_ = t_n_doc;
    n_nodes_ = node_starts_.size();
    block_size_ = t_block_size;
    // A policy that does not read the storing factor builds the same core for
    // every value of it; record none, so the cores are the same bytes too.
    storing_factor_ = PolicyReadsStoringFactor(t_policy) ? t_storing_factor : 0.0f;
    policy_ = t_policy;
    rebindRankSelect();
  }

  // Serialized in component order (see CacheComponents below): the tree, which
  // depends on the block size only; the nodes the policy selects; their lists.
  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    std::size_t bytes = 0;
    bytes += node_starts_.serialize(out, child, "node_starts");
    bytes += node_ends_.serialize(out, child, "node_ends");
    bytes += first_child_.serialize(out, child, "first_child");
    bytes += next_sibling_.serialize(out, child, "next_sibling");
    bytes += sdsl::write_member(n_nodes_, out, child, "n_nodes");
    bytes += sdsl::write_member(block_size_, out, child, "block_size");
    bytes += sdsl::write_member(n_doc_, out, child, "n_doc");
    bytes += selected_marker_.serialize(out, child, "selected_marker");
    auto policy_id = static_cast<uint8_t>(policy_);
    bytes += sdsl::write_member(policy_id, out, child, "policy");
    bytes += sdsl::write_member(storing_factor_, out, child, "storing_factor");
    bytes += stored_sets_.serialize(out, child, "stored_sets");
    sdsl::structure_tree::add_size(child, bytes);
    return bytes;
  }

  void load(std::istream& in) {
    node_starts_.load(in);
    node_ends_.load(in);
    first_child_.load(in);
    next_sibling_.load(in);
    sdsl::read_member(n_nodes_, in);
    sdsl::read_member(block_size_, in);
    sdsl::read_member(n_doc_, in);
    selected_marker_.load(in);
    uint8_t policy_id = 0;
    sdsl::read_member(policy_id, in);
    policy_ = static_cast<StoragePolicy>(policy_id);
    sdsl::read_member(storing_factor_, in);
    stored_sets_.load(in);
    rebindRankSelect();
  }

  // The format cached before 2026-10 (marker first, scalars last, the storing
  // factor recorded under every policy). Read only by the cache migration tool,
  // which then normalizes the storing factor as Assemble does.
  void loadPre202610(std::istream& in) {
    selected_marker_.load(in);
    node_starts_.load(in);
    node_ends_.load(in);
    first_child_.load(in);
    next_sibling_.load(in);
    stored_sets_.load(in);
    sdsl::read_member(n_doc_, in);
    sdsl::read_member(n_nodes_, in);
    sdsl::read_member(block_size_, in);
    sdsl::read_member(storing_factor_, in);
    uint8_t policy_id = 0;
    sdsl::read_member(policy_id, in);
    policy_ = static_cast<StoragePolicy>(policy_id);
    if (!PolicyReadsStoringFactor(policy_)) storing_factor_ = 0.0f;
    rebindRankSelect();
  }

  // Bytes of the tree and selection components (the lists are the rest).
  std::size_t TreeComponentSize() const {
    return sdsl::size_in_bytes(node_starts_) + sdsl::size_in_bytes(node_ends_) + sdsl::size_in_bytes(first_child_) +
           sdsl::size_in_bytes(next_sibling_) + sizeof(n_nodes_) + sizeof(block_size_) + sizeof(n_doc_);
  }
  std::size_t SelectionComponentSize() const {
    return sdsl::size_in_bytes(selected_marker_) + sizeof(uint8_t) + sizeof(storing_factor_);
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    append(r, "selected_marker", sdsl::size_in_bytes(selected_marker_));
    append(r, "node_starts", sdsl::size_in_bytes(node_starts_));
    append(r, "node_ends", sdsl::size_in_bytes(node_ends_));
    append(r, "first_child", sdsl::size_in_bytes(first_child_));
    append(r, "next_sibling", sdsl::size_in_bytes(next_sibling_));
    auto stored = stored_sets_.GetSizeReport();
    for (const auto& f : stored) append(r, "stored_sets_" + f.name, f.bytes);
    return r;
  }

  std::size_t n_doc() const { return n_doc_; }
  std::size_t n_nodes() const { return n_nodes_; }
  uint32_t block_size() const { return block_size_; }
  float storing_factor() const { return storing_factor_; }
  StoragePolicy policy() const { return policy_; }

 private:
  // Sentinel for "no child / no sibling". Equals n_nodes_ so
  // first_child_[id] / next_sibling_[id] are sized to fit n_nodes_.
  std::size_t kSentinel() const { return n_nodes_; }

  // Re-point rank/select supports at selected_marker_ after any move/load.
  // Per the project's loadItemPtr gotcha: SDSL rank/select hold raw
  // pointers into the underlying bitvector and become dangling after
  // moves. Re-binding from the now-stable member address is safe.
  void rebindRankSelect() {
    selected_rank_ = TBvRank(&selected_marker_);
    selected_select_ = TBvSelect(&selected_marker_);
  }

  TBitvector selected_marker_{};
  TBvRank selected_rank_{};
  TBvSelect selected_select_{};
  TIntVector node_starts_{};
  TIntVector node_ends_{};
  TIntVector first_child_{};
  TIntVector next_sibling_{};
  TStoredSetCodec stored_sets_{};
  std::size_t n_doc_ = 0;
  std::size_t n_nodes_ = 0;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4.0f;
  StoragePolicy policy_ = StoragePolicy::OccurrenceWeighted;
};

// TCore with its stored-set codec replaced by TCodec.
template <typename TCore, typename TCodec>
struct WithCodec;
template <typename TBitvector, typename TBvRank, typename TBvSelect, typename TIntVector, typename TOld, typename TCodec>
struct WithCodec<PDLTreeCore<TBitvector, TBvRank, TBvSelect, TIntVector, TOld>, TCodec> {
  using type = PDLTreeCore<TBitvector, TBvRank, TBvSelect, TIntVector, TCodec>;
};

// Prefix of the PDL components that depend on the selection: the block size,
// the policy and, for a policy that reads it, the storing factor.
inline std::string PdlSelectionPrefix(const JSON& t_keys,
                                      uint32_t t_block_size,
                                      float t_storing_factor,
                                      StoragePolicy t_policy) {
  const auto& terms = t_keys[conf::kPolicy];
  switch (t_policy) {
    case StoragePolicy::OccurrenceWeighted:
      return PrefixedKey(t_keys, conf::kBlkSfPol, "", t_block_size, t_storing_factor,
                         terms[conf::kOccW].get<std::string>());
    case StoragePolicy::ListWeighted:
      return PrefixedKey(t_keys, conf::kBlkSfPol, "", t_block_size, t_storing_factor,
                         terms[conf::kListW].get<std::string>());
    case StoragePolicy::LeavesOnly:
      return PrefixedKey(t_keys, conf::kBlkPol, "", t_block_size, terms[conf::kLeaves].get<std::string>());
    case StoragePolicy::StoreAllInternal:
      return PrefixedKey(t_keys, conf::kBlkPol, "", t_block_size, terms[conf::kAllNodes].get<std::string>());
  }
  throw std::invalid_argument("unknown PDL storage policy");
}

inline std::string PdlSelectionKey(const JSON& t_keys,
                                   std::string_view t_name,
                                   uint32_t t_block_size,
                                   float t_storing_factor,
                                   StoragePolicy t_policy) {
  return PdlSelectionPrefix(t_keys, t_block_size, t_storing_factor, t_policy) + t_keys[t_name].get<std::string>();
}

// Cache components of a PDL core, in the order it serializes them: the tree
// (block size only, shared by every policy, storing factor and codec), the
// selected nodes (shared by the codecs) and the lists in the core's codec.
template <typename TBitvector, typename TBvRank, typename TBvSelect, typename TIntVector, typename TStoredSetCodec>
std::vector<Component> CacheComponents(
    const PDLTreeCore<TBitvector, TBvRank, TBvSelect, TIntVector, TStoredSetCodec>& t_core,
    const JSON& t_keys,
    uint32_t t_block_size,
    float t_storing_factor,
    StoragePolicy t_policy) {
  const auto tree = t_core.TreeComponentSize();
  const auto selection = t_core.SelectionComponentSize();
  return {
      {PrefixedKey(t_keys, conf::kBlk, t_keys[conf::kPdlTree].get<std::string>(), t_block_size),
       TypeHash<TIntVector>(), tree},
      {PdlSelectionKey(t_keys, conf::kPdlSelection, t_block_size, t_storing_factor, t_policy),
       TypeHash<TBitvector>(), selection},
      {PdlSelectionKey(t_keys, TStoredSetCodec::kListsKey, t_block_size, t_storing_factor, t_policy),
       TypeHash<TStoredSetCodec>(), SerializedSize(t_core) - tree - selection},
  };
}

template <typename TBitvector, typename TBvRank, typename TBvSelect,
          typename TIntVector, typename TStoredSetCodec>
void collectSizes(SizeReport& out,
                  const PDLTreeCore<TBitvector, TBvRank, TBvSelect,
                                    TIntVector, TStoredSetCodec>& core,
                  const std::string& prefix = "") {
  for (const auto& f : core.GetSizeReport()) append(out, prefix + f.name, f.bytes);
}

}  // namespace dret::pdl
