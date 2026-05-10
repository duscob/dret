//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// Skeleton for the PDL sparse-suffix-tree core. Construction logic, real
// cover computation, and codec wiring land in later tasks (see
// docs/pdl_indexes_tasks.md, Tasks 6–13). For now this header just compiles
// with default template arguments so dependent classes (Tasks 23–25) can be
// declared in parallel.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "../size_report.h"

namespace dret::pdl {

enum class StoragePolicy : uint8_t {
  OriginalDrl = 0,
  StoreAllInternal = 1,
  LeavesOnly = 2,
};

// Placeholder codec satisfying the minimum API used by PDLTreeCore. Replaced
// by Plain / RP / BC codecs in Tasks 15–17.
struct NullCodec {
  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    sdsl::structure_tree::add_size(child, 0);
    return 0;
  }
  void load(std::istream& /*in*/) {}
};

template <typename TBitvector = sdsl::sd_vector<>,
          typename TBvRank = typename TBitvector::rank_1_type,
          typename TBvSelect = typename TBitvector::select_1_type,
          typename TIntVector = sdsl::int_vector<>,
          typename TStoredSetCodec = NullCodec>
class PDLTreeCore {
 public:
  using TDocId = std::size_t;

  PDLTreeCore() = default;

  // Single-range cover; matches DLSampledTreeScheme::computeCover. Real impl
  // in Task 12; the skeleton returns an empty selected-node set so callers
  // fall back to raw-range retrieval over the whole [sp, ep).
  std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> computeCover(
      std::size_t t_sp,
      std::size_t /*t_ep*/) const {
    return {{t_sp, t_sp}, {}};
  }

  // Multi-range cover; needed because non-OriginalDrl storage policies can
  // leave non-contiguous gaps in [sp, ep). The default delegates to
  // computeCover and is safe for any single-range subclass; PDLTreeCore will
  // override this once Task 12 wires real navigation in.
  void computeCoverFull(std::size_t t_sp,
                        std::size_t t_ep,
                        std::vector<std::pair<std::size_t, std::size_t>>& t_raw_ranges,
                        std::vector<std::size_t>& t_nodes) const {
    auto [range, nodes] = computeCover(t_sp, t_ep);
    if (nodes.empty()) {
      t_raw_ranges.emplace_back(t_sp, t_ep);
    } else {
      if (t_sp < range.first) t_raw_ranges.emplace_back(t_sp, range.first);
      if (range.second < t_ep) t_raw_ranges.emplace_back(range.second, t_ep);
    }
    t_nodes = std::move(nodes);
  }

  std::vector<TDocId> getDocSet(std::size_t /*t_node_id*/) const {
    return {};
  }

  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    std::size_t bytes = 0;
    bytes += selected_marker_.serialize(out, child, "selected_marker");
    bytes += node_starts_.serialize(out, child, "node_starts");
    bytes += node_ends_.serialize(out, child, "node_ends");
    bytes += stored_sets_.serialize(out, child, "stored_sets");
    bytes += sdsl::write_member(n_doc_, out, child, "n_doc");
    bytes += sdsl::write_member(block_size_, out, child, "block_size");
    bytes += sdsl::write_member(storing_factor_, out, child, "storing_factor");
    auto policy_id = static_cast<uint8_t>(policy_);
    bytes += sdsl::write_member(policy_id, out, child, "policy");
    sdsl::structure_tree::add_size(child, bytes);
    return bytes;
  }

  void load(std::istream& in) {
    selected_marker_.load(in);
    node_starts_.load(in);
    node_ends_.load(in);
    stored_sets_.load(in);
    sdsl::read_member(n_doc_, in);
    sdsl::read_member(block_size_, in);
    sdsl::read_member(storing_factor_, in);
    uint8_t policy_id = 0;
    sdsl::read_member(policy_id, in);
    policy_ = static_cast<StoragePolicy>(policy_id);
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    append(r, "selected_marker", sdsl::size_in_bytes(selected_marker_));
    append(r, "node_starts", sdsl::size_in_bytes(node_starts_));
    append(r, "node_ends", sdsl::size_in_bytes(node_ends_));
    return r;
  }

  std::size_t n_doc() const { return n_doc_; }
  uint32_t block_size() const { return block_size_; }
  float storing_factor() const { return storing_factor_; }
  StoragePolicy policy() const { return policy_; }

 private:
  TBitvector selected_marker_{};
  TBvRank selected_rank_{};
  TBvSelect selected_select_{};
  TIntVector node_starts_{};
  TIntVector node_ends_{};
  TStoredSetCodec stored_sets_{};
  std::size_t n_doc_ = 0;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4.0f;
  StoragePolicy policy_ = StoragePolicy::OriginalDrl;
};

}  // namespace dret::pdl
