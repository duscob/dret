//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// Common API for PDL stored-set codecs (Tasks 14-17). Every codec
// stores one document set per selected node, indexed by codec slot
// (the rank-1 of the node's id in PDLTreeCore::selected_marker_).
//
// Concrete codecs land in later tasks:
//   - PlainCodec  (Task 15) — sorted doc-id payload per slot.
//   - RPCodec     (Task 16) — RePair-compressed payload, recursively expanded.
//   - BCCodec     (Task 17) — dictionary-coded blocks; design fresh.
//
// All three carry the all-doc sentinel: at build time a node whose
// docs.size() reached n_doc is recorded as contains_all=true with no
// payload; at expansion time that slot yields 0..n_doc-1 via
// ExpandAllDoc.
//

#pragma once

#include <concepts>
#include <cstddef>
#include <functional>
#include <istream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <grammar/slp_metadata.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include "../size_report.h"

namespace dret::pdl {

// One stored-set entry as seen by the codec at build time. contains_all
// = true means the all-doc sentinel: the codec stores nothing for the
// slot and Expand emits 0..n_doc-1. Otherwise docs is sorted-unique
// zero-based ids; the codec is free to compress them however it likes.
struct StoredSet {
  bool contains_all = false;
  std::vector<std::size_t> docs;
};

// Helper: emit 0..n_doc-1 into a callback. Used by every codec on the
// all-doc sentinel path.
template <typename TReport>
void ExpandAllDoc(std::size_t t_n_doc, TReport&& t_report) {
  for (std::size_t i = 0; i < t_n_doc; ++i) t_report(i);
}

// Codec interface. The Build method's signature varies per codec
// (different internal containers), so it's documented here rather than
// concept-enforced; Tasks 15-17 each accept a (n_slots, get_set_at)
// pair, where get_set_at(slot) returns a StoredSet by const reference
// or by value.
//
// The concept covers the query-time contract (Expand) and the
// persistence contract (serialize / load). It uses a concrete std::
// function as the callback in the concept body so the requires-clause
// stays checkable; codecs are free to implement Expand as a template
// method that accepts any callable for performance.
template <typename T>
concept SetCodec = std::default_initializable<T>
    && requires(T& mut, const T& cnst, std::ostream& os, std::istream& is) {
         { cnst.Expand(std::size_t{}, std::size_t{},
                       std::function<void(std::size_t)>{}) } -> std::same_as<void>;
         { cnst.serialize(os) } -> std::convertible_to<std::size_t>;
         { mut.load(is) } -> std::same_as<void>;
       };

// Reference codec exemplifying the SetCodec interface. Stores nothing;
// Build is a no-op; Expand always yields the all-doc set. Useful as a
// compile-time sanity check for the concept and as a no-op stand-in in
// tests that don't care about real codec output.
class DummyCodec {
 public:
  template <typename TGetSetAt>
  void Build(std::size_t /*t_n_slots*/, TGetSetAt&& /*t_get_set_at*/) {}

  template <typename TReport>
  void Expand(std::size_t /*t_slot*/, std::size_t t_n_doc, TReport&& t_report) const {
    ExpandAllDoc(t_n_doc, std::forward<TReport>(t_report));
  }

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

static_assert(SetCodec<DummyCodec>,
              "DummyCodec must satisfy the SetCodec concept (Task 14 acceptance)");

// PlainCodec — sorted-doc-id payload per slot, backed by grammar::Chunks
// (a flat objects vector + per-chunk start positions). Templated on the
// object and position container types so Task 18 can swap them for
// compressed alternatives.
//
// All-doc sentinel: stored as a single-element chunk holding the value
// n_doc (one past the maximum legal doc id). At Expand time a chunk that
// is exactly [n_doc] dispatches to ExpandAllDoc; everything else
// iterates the chunk verbatim. The sentinel is recognised at expansion,
// not at construction time, so the codec doesn't need to store n_doc
// itself — the caller passes it via the Expand argument.
//
// Build signature: (n_slots, get_set_at, n_doc). get_set_at(slot) must
// return a StoredSet with sorted-unique zero-based docs.
template <typename TObjContainer = sdsl::int_vector<>,
          typename TPosContainer = sdsl::int_vector<>>
class PlainCodec {
 public:
  using TStoredChunks = grammar::Chunks<TObjContainer, TPosContainer>;

  PlainCodec() = default;

  template <typename TGetSetAt>
  void Build(std::size_t t_n_slots, TGetSetAt&& t_get_set_at, std::size_t t_n_doc) {
    // Build into a std::vector-backed Chunks (which has push_back), then
    // copy-construct into the (possibly compressed) TObjContainer/TPosContainer.
    grammar::Chunks<std::vector<std::size_t>, std::vector<std::size_t>> tmp;
    for (std::size_t s = 0; s < t_n_slots; ++s) {
      StoredSet set = t_get_set_at(s);
      if (set.contains_all) {
        tmp.Insert(t_n_doc);
      } else {
        tmp.Insert(set.docs.begin(), set.docs.end());
      }
    }
    chunks_ = TStoredChunks(tmp);
  }

  template <typename TReport>
  void Expand(std::size_t t_slot, std::size_t t_n_doc, TReport&& t_report) const {
    auto chunk = chunks_[t_slot + 1];  // grammar::Chunks is 1-indexed.
    if (chunk.size() == 1
        && static_cast<std::size_t>(*chunk.begin()) == t_n_doc) {
      ExpandAllDoc(t_n_doc, std::forward<TReport>(t_report));
      return;
    }
    for (auto it = chunk.begin(); it != chunk.end(); ++it) {
      t_report(static_cast<std::size_t>(*it));
    }
  }

  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    auto bytes = chunks_.serialize(out, child, "chunks");
    sdsl::structure_tree::add_size(child, bytes);
    return bytes;
  }

  void load(std::istream& in) {
    chunks_.load(in);
  }

  std::size_t n_slots() const { return chunks_.size(); }

  SizeReport GetSizeReport() const {
    SizeReport r;
    append(r, "chunks_objs", sdsl::size_in_bytes(chunks_.GetObjects()));
    append(r, "chunks_pos", sdsl::size_in_bytes(chunks_.GetChunksPositions()));
    return r;
  }

 private:
  TStoredChunks chunks_;
};

static_assert(SetCodec<PlainCodec<>>,
              "PlainCodec<> must satisfy the SetCodec concept");

}  // namespace dret::pdl
