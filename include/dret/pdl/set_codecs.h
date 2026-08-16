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

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <istream>
#include <iterator>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_metadata.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "dret/config.h"
#include "dret/size_report.h"

#ifndef VNMEXTRACT_EXE
#define VNMEXTRACT_EXE "vnmextract"
#endif

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
  static constexpr bool kExpandsSorted = true;

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
  static constexpr std::string_view kVariantKey = conf::kPlain;

  // Expand emits doc ids in ascending order, so PDLTreeCore::getDocSet can
  // hand the result straight to std::set_union in the merge functor.
  static constexpr bool kExpandsSorted = true;

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
    // Bit-compress the int_vector-backed objs/pos on convert (no-op for other
    // containers). Without an action the chunks stay at the default 64-bit width.
    auto bit_compress = [](auto& v) {
      if constexpr (std::is_same_v<std::decay_t<decltype(v)>, sdsl::int_vector<>>)
        sdsl::util::bit_compress(v);
    };
    chunks_ = TStoredChunks(tmp, bit_compress, bit_compress);
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

// RPCodec — RePair-compressed stored sets via grammar::GCChunks. All
// slots share a single RePair grammar (built over the concatenated
// per-slot doc-id sequences) plus per-slot compressed compact sequences;
// a query expands its slot's compact sequence by recursively walking
// the shared SLP. Ports the recursive-expansion idea from
// drl/src/pdlrp.cpp:202-236, but reuses grammar::RePairEncoder<false>
// (the encoder GCDA already uses for its chunk grammar) instead of
// reinventing RePair.
//
// All-doc sentinel handling is identical to PlainCodec: at Build time
// the slot's "set" is a single n_doc; at Expand time the expanded
// chunk vector with a single n_doc dispatches to ExpandAllDoc. RePair
// might compress runs of n_doc into rules, but each slot's chunk
// expansion still yields a single n_doc per all-doc slot.
//
// Templated on the underlying SLP and per-slot Chunks types so Task 18
// can swap them; defaults match GCDA's chunk grammar setup so the same
// type-hashed cache entries can be shared if needed.
// TSLP defaults to grammar::SLP with sdsl::int_vector<> rule + leaf containers
// (bit-compressed) rather than grammar::SLP<>'s default std::vector — same
// precedent as the GCDA combined/differential base-grammar fix
// (docs/container_audit.md). The Stage-3 generic-action build path
// (compress_if_iv lambda below) calls sdsl::util::bit_compress on
// sdsl::int_vector<> fields, so this default is the bit-compressed stored
// form. (Existing rp cache files keyed by the prior std::vector type-hash
// will be rebuilt on first run.)
template <typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
          typename TChunks = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>
class RPCodec {
 public:
  static constexpr std::string_view kVariantKey = conf::kRP;

  // Expand emits doc ids in ascending order, so PDLTreeCore::getDocSet can
  // hand the result straight to std::set_union in the merge functor.
  static constexpr bool kExpandsSorted = true;

  using TStoredChunks = grammar::GCChunks<TSLP, /*kExpand=*/true, TChunks>;

  RPCodec() = default;

  template <typename TGetSetAt>
  void Build(std::size_t t_n_slots, TGetSetAt&& t_get_set_at, std::size_t t_n_doc) {
    // Stage 1: build a temporary uncompressed Chunks holding all slots.
    grammar::Chunks<std::vector<std::size_t>, std::vector<std::size_t>> tmp;
    for (std::size_t s = 0; s < t_n_slots; ++s) {
      StoredSet set = t_get_set_at(s);
      if (set.contains_all) {
        tmp.Insert(t_n_doc);
      } else {
        tmp.Insert(set.docs.begin(), set.docs.end());
      }
    }

    // Stage 2: encode via RePair into a default (std::vector-backed)
    // intermediate GCChunks. The detour through std::vector storage is
    // necessary because grammar/slp_metadata.h:202 does
    // back_inserter(objs_) unqualified — ADL finds std::back_inserter only
    // when objs_ lives in std (i.e., std::vector). Mirrors the staged
    // construct() path GCDA uses (doc_list/doc_list_gcda.h:474-516).
    grammar::GCChunks<TSLP> intermediate;
    grammar::RePairEncoder<false> encoder;
    const auto& objs = tmp.GetObjects();
    intermediate.Compute(objs.begin(), objs.end(), tmp, encoder);

    // Stage 3: copy-convert into the configured TStoredChunks. The
    // generic action bit-compresses sdsl::int_vector<> fields when the
    // configured TSLP/TChunks use them, and is a no-op for std::vector
    // fields (e.g., grammar::SLP<>'s default storage).
    auto compress_if_iv = [](auto& v) {
      if constexpr (std::is_same_v<std::remove_reference_t<decltype(v)>,
                                   sdsl::int_vector<>>) {
        sdsl::util::bit_compress(v);
      }
    };
    chunks_ = TStoredChunks(intermediate, compress_if_iv, compress_if_iv,
                            compress_if_iv, compress_if_iv);
  }

  template <typename TReport>
  void Expand(std::size_t t_slot, std::size_t t_n_doc, TReport&& t_report) const {
    auto v = chunks_[t_slot + 1];  // expanded set; grammar::Chunks is 1-indexed.
    if (v.size() == 1 && static_cast<std::size_t>(v[0]) == t_n_doc) {
      ExpandAllDoc(t_n_doc, std::forward<TReport>(t_report));
      return;
    }
    for (auto x : v) t_report(static_cast<std::size_t>(x));
  }

  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    auto bytes = chunks_.serialize(out, child, "gc_chunks");
    sdsl::structure_tree::add_size(child, bytes);
    return bytes;
  }

  void load(std::istream& in) { chunks_.load(in); }

  std::size_t n_slots() const { return chunks_.size(); }

  SizeReport GetSizeReport() const {
    SizeReport r;
    append(r, "gc_chunks", sdsl::size_in_bytes(chunks_));
    return r;
  }

 private:
  TStoredChunks chunks_;
};

static_assert(SetCodec<RPCodec<>>,
              "RPCodec<> must satisfy the SetCodec concept");

// BCCodec — Block / biClique-coded stored sets, mirroring drl's PDL-BC
// (CSA::DocArray, drl/external/rlcsa/docarray.{h,cpp}). Each slot's set
// is stored as a "block": a flat sequence of values where
//
//   value < n_doc                  → terminal doc id (0-based)
//   n_doc <= value < n_doc + R     → reference to rule (value - n_doc)
//   value >= n_doc + R             → all-doc sentinel
//
// where R is the number of rules. Rules are flat doc-id sequences shared
// across blocks; they are produced by Cecilia Hernandez's web-graph
// dense-subgraph extractor (vnmextract, see external/dsextract/) on the
// bipartite (slot, doc) graph. drl's design assumes bicliques fully
// cover every block; we relax that by allowing blocks to mix rule-refs
// with verbatim doc-ids so small fixtures (where vnmextract finds no
// bicliques) still encode correctly. This stays compatible with drl's
// expansion algorithm — the dispatch above doesn't care about the order
// of doc-ids and rule-refs inside a block.
//
// Build invokes vnmextract via std::system using the VNMEXTRACT_EXE
// compile-define (set by ConfigDSExtract.cmake). It writes a binary
// .graph file in drl's document_graph.cpp format (header: u32 numNodes,
// u32 numEdges; body: int32 stream where negative = -node_id starts a
// slot, positive = doc+1 is an outedge), then parses
// `<prefix>-biclique-it-N.txt` files (text: `block_ids - doc_ids` per
// row) into rules + per-slot blocks.
template <typename TBitvector = sdsl::sd_vector<>,
          typename TIntVector = sdsl::int_vector<>>
class BCCodec {
 public:
  static constexpr std::string_view kVariantKey = conf::kBC;

  static constexpr bool kExpandsSorted = false;

  using TSelect1 = typename TBitvector::select_1_type;

  // vnmextract mining parameters, matching its positional argv:
  //
  //   vnmextract <graph> <format> <shingle_size> <iters> <bcsizes> <out> <num_hashes>
  //
  // format=1 (binary) and shingle_size=1 are fixed by runVnmextract; the other
  // three are these. Names follow vnmextract's own usage string rather than
  // what they might sound like:
  //
  //   min_bicliques (argv[4], `iters` upstream) is NOT a node-degree cutoff.
  //     It is a minimum-yield stopping rule. Each pass mines at the current
  //     bcsize, then compares the number of bicliques found against it; on a
  //     smaller yield the miner drops to the next bcsize, and stops once the
  //     list is exhausted. See the loop at the end of vnmextract.cpp main().
  //   bcsizes (argv[5]) is the descending list of target biclique sizes.
  //   num_hashes (argv[7]) is the number of min-hash functions used for
  //     shingling, not a count of shingles.
  //
  // Defaults are drl's production values (drl/external/doclist-env/build_pdl):
  // 500, 5000,500,100,50,30,15, and 4. drl used 1000 and 8 for its single
  // largest (1.09 GB) collection, enwiki-big.
  //
  // These were previously hardcoded to 1 / 10,5,2 so that small test fixtures
  // would still yield bicliques — but nothing confined that to fixtures, and it
  // applied to real collections too. min_bicliques=1 means every pass keeps
  // mining at the current size until one yields literally zero, and mining down
  // to bcsize=2 continues long past where drl stopped. Measured: 659 mining
  // iterations on a 100 MB collection, 17.6x slower on `page` for 0.15% less
  // space, and >55 h without finishing one DA row on a 671 MB collection that
  // takes 218 s at these defaults. See docs/bug_pdlbc_parameters.md.
  //
  // Fixtures that need bicliques out of a tiny input must opt in explicitly via
  // SetMiningParams — see PDLBCCodec.RulesAndBlocksRoundTripWhenBicliquesExtracted.
  struct MiningParams {
    std::string min_bicliques = "500";
    std::string bcsizes = "5000,500,100,50,30,15";
    std::string num_hashes = "4";
  };

  BCCodec() = default;

  void SetMiningParams(MiningParams t_params) { mining_ = std::move(t_params); }

  const MiningParams& mining_params() const { return mining_; }

  template <typename TGetSetAt>
  void Build(std::size_t t_n_slots, TGetSetAt&& t_get_set_at, std::size_t t_n_doc) {
    n_doc_ = t_n_doc;
    n_slots_ = t_n_slots;
    n_rules_ = 0;
    blocks_ = TIntVector{};
    rules_ = TIntVector{};
    block_borders_ = TBitvector{};
    rule_borders_ = TBitvector{};
    if (t_n_slots == 0) {
      rebindSelect();
      return;
    }

    // Materialise the per-slot StoredSets so we can pass over them twice
    // (once to write the .graph, once to compute remaining-uncovered
    // docs after biclique extraction).
    std::vector<StoredSet> slot_sets;
    slot_sets.reserve(t_n_slots);
    for (std::size_t s = 0; s < t_n_slots; ++s) slot_sets.push_back(t_get_set_at(s));

    // Stage 1: write the .graph binary file. Skip contains_all slots and
    // empty slots — vnmextract takes only edges, and we represent the
    // sentinel / empty cases directly in the block content later.
    auto tmpdir = make_temp_dir();
    auto graph_path = (tmpdir / "input.graph").string();
    auto out_prefix = (tmpdir / "out").string();

    std::vector<std::pair<std::int32_t, std::vector<std::uint32_t>>> graph_nodes;
    graph_nodes.reserve(t_n_slots);
    std::uint32_t total_edges = 0;
    for (std::size_t s = 0; s < t_n_slots; ++s) {
      const auto& set = slot_sets[s];
      // Singletons are excluded from the graph, matching the reference
      // implementation (rlcsa document_graph.cpp:209, which writes them to a
      // separate .singletons file with the note "The presence of singletons
      // greatly slows down the search for bicliques"). A degree-1 node cannot
      // participate in a biclique -- that needs at least two documents on the
      // document side -- so it is pure search-space dilution.
      //
      // Measured on version_0100_100 at (1024,8): 618,288 of 2,428,115 graph
      // nodes (25.5%) were singletons. They inflate the row count that
      // vnmextract shingles over, so min-hash clusters fill up with nodes that
      // can never contribute a biclique -- which both slows every pass and
      // lowers its yield. See docs/bug_pdlbc_parameters.md.
      //
      // Stage 3 below is unaffected: a slot absent from the graph simply has
      // no rule refs and emits its documents verbatim, which is exactly how
      // the reference treated its .singletons entries.
      if (set.contains_all || set.docs.empty() || set.docs.size() == 1) continue;
      // drl convention: slot node id = n_doc + 1 + slot, written negative;
      // doc id is 1-based (doc + 1).
      auto node_id = static_cast<std::int32_t>(t_n_doc + 1 + s);
      std::vector<std::uint32_t> edges;
      edges.reserve(set.docs.size());
      for (auto d : set.docs) edges.push_back(static_cast<std::uint32_t>(d + 1));
      total_edges += edges.size();
      graph_nodes.emplace_back(-node_id, std::move(edges));
    }

    if (!graph_nodes.empty()) {
      writeGraphFile(graph_path, graph_nodes, total_edges);
      runVnmextract(graph_path, out_prefix);
    }

    // Stage 2: parse all biclique files in the temp dir.
    std::vector<std::vector<std::size_t>> rule_docs;        // per-rule doc ids
    std::vector<std::vector<std::size_t>> slot_rule_refs;   // per-slot rule indices
    slot_rule_refs.resize(t_n_slots);
    parseBicliqueFiles(out_prefix, t_n_doc, t_n_slots, rule_docs, slot_rule_refs);

    // Stage 3: for each slot, pick a non-overlapping subset of its
    // rule-refs (greedy: take rules in the order they appear; drop a
    // rule if it would duplicate any doc already covered) and emit
    // remaining-uncovered docs verbatim. This keeps Expand correct
    // regardless of how aggressive vnmextract was — and also matches
    // PlainCodec output exactly when no rules apply.
    std::vector<std::vector<std::size_t>> slot_blocks(t_n_slots);
    for (std::size_t s = 0; s < t_n_slots; ++s) {
      const auto& set = slot_sets[s];
      auto& block = slot_blocks[s];
      if (set.contains_all) {
        block.push_back(t_n_doc + rule_docs.size());  // sentinel = maxInteger
        continue;
      }
      std::unordered_set<std::size_t> covered;
      std::vector<std::size_t> kept_rules;
      for (auto rid : slot_rule_refs[s]) {
        const auto& rdocs = rule_docs[rid];
        bool overlaps = false;
        for (auto d : rdocs) {
          if (covered.count(d)) { overlaps = true; break; }
        }
        if (overlaps) continue;
        for (auto d : rdocs) covered.insert(d);
        kept_rules.push_back(rid);
      }
      // Emit the kept rule refs first, then verbatim docs not covered by
      // any kept rule. Order doesn't matter for correctness — the caller
      // dedupes / sorts.
      for (auto rid : kept_rules) block.push_back(t_n_doc + rid);
      for (auto d : set.docs) {
        if (!covered.count(d)) block.push_back(d);
      }
    }

    // Stage 4: pack rules and blocks into TIntVector + sd_vector borders.
    n_rules_ = rule_docs.size();
    packRules(rule_docs, t_n_doc);
    packBlocks(slot_blocks);

    cleanupTempDir(tmpdir);
    rebindSelect();
  }

  template <typename TReport>
  void Expand(std::size_t t_slot, std::size_t t_n_doc, TReport&& t_report) const {
    if (t_slot >= n_slots_) return;
    const std::size_t b_from = block_select_(t_slot + 1);
    const std::size_t b_to =
        (t_slot + 1 < n_slots_) ? block_select_(t_slot + 2) : blocks_.size();
    const std::size_t max_int = t_n_doc + n_rules_;
    for (std::size_t i = b_from; i < b_to; ++i) {
      const auto val = static_cast<std::size_t>(blocks_[i]);
      if (val >= max_int) {
        ExpandAllDoc(t_n_doc, std::forward<TReport>(t_report));
        return;
      }
      if (val < t_n_doc) {
        t_report(val);
        continue;
      }
      // Rule reference.
      const std::size_t rule_idx = val - t_n_doc;
      const std::size_t r_from = rule_select_(rule_idx + 1);
      const std::size_t r_to = (rule_idx + 1 < n_rules_) ? rule_select_(rule_idx + 2)
                                                         : rules_.size();
      for (std::size_t j = r_from; j < r_to; ++j) {
        const auto rval = static_cast<std::size_t>(rules_[j]);
        if (rval >= t_n_doc) {
          ExpandAllDoc(t_n_doc, std::forward<TReport>(t_report));
          return;
        }
        t_report(rval);
      }
    }
  }

  std::size_t n_slots() const { return n_slots_; }
  std::size_t n_rules() const { return n_rules_; }
  std::size_t n_doc() const { return n_doc_; }

  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    std::size_t bytes = 0;
    bytes += sdsl::write_member(n_doc_, out, child, "n_doc");
    bytes += sdsl::write_member(n_slots_, out, child, "n_slots");
    bytes += sdsl::write_member(n_rules_, out, child, "n_rules");
    bytes += blocks_.serialize(out, child, "blocks");
    bytes += block_borders_.serialize(out, child, "block_borders");
    bytes += rules_.serialize(out, child, "rules");
    bytes += rule_borders_.serialize(out, child, "rule_borders");
    sdsl::structure_tree::add_size(child, bytes);
    return bytes;
  }

  void load(std::istream& in) {
    sdsl::read_member(n_doc_, in);
    sdsl::read_member(n_slots_, in);
    sdsl::read_member(n_rules_, in);
    blocks_.load(in);
    block_borders_.load(in);
    rules_.load(in);
    rule_borders_.load(in);
    rebindSelect();
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    append(r, "blocks", sdsl::size_in_bytes(blocks_));
    append(r, "block_borders", sdsl::size_in_bytes(block_borders_));
    append(r, "rules", sdsl::size_in_bytes(rules_));
    append(r, "rule_borders", sdsl::size_in_bytes(rule_borders_));
    return r;
  }

 private:
  // Re-point select supports at the loaded/assembled bitvectors. Same
  // dangling-pointer guard as PDLTreeCore::rebindRankSelect.
  void rebindSelect() {
    block_select_ = TSelect1(&block_borders_);
    rule_select_ = TSelect1(&rule_borders_);
  }

  static std::filesystem::path make_temp_dir() {
    auto base = std::filesystem::temp_directory_path() / "dret_bc_XXXXXX";
    std::string tmpl = base.string();
    if (mkdtemp(tmpl.data()) == nullptr) {
      throw std::runtime_error("BCCodec: mkdtemp failed");
    }
    return std::filesystem::path(tmpl);
  }

  static void cleanupTempDir(const std::filesystem::path& dir) {
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
  }

  static void writeGraphFile(
      const std::string& path,
      const std::vector<std::pair<std::int32_t, std::vector<std::uint32_t>>>& nodes,
      std::uint32_t total_edges) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("BCCodec: cannot open .graph for write: " + path);
    auto num_nodes = static_cast<std::uint32_t>(nodes.size());
    out.write(reinterpret_cast<const char*>(&num_nodes), sizeof(num_nodes));
    out.write(reinterpret_cast<const char*>(&total_edges), sizeof(total_edges));
    for (const auto& [neg_node_id, edges] : nodes) {
      out.write(reinterpret_cast<const char*>(&neg_node_id), sizeof(neg_node_id));
      for (auto e : edges) out.write(reinterpret_cast<const char*>(&e), sizeof(e));
    }
  }

  // Runs the external biclique miner. Throws if it cannot be run at all.
  //
  // This used to be `(void)rc` with stdout and stderr sent to /dev/null, on the
  // reasoning that vnmextract legitimately emits no biclique files when there
  // is nothing compressible. That is true, but it made three very different
  // situations indistinguishable: no bicliques found, the miner crashing, and
  // the binary not existing at all. The last one is not hypothetical --
  // `bm_doc_list` had no build dependency on the dsextract external project, so
  // any fresh build tree produced a binary whose VNMEXTRACT_EXE pointed at a
  // missing file. Every PDL-BC build then silently reported zero rules, which
  // is indistinguishable from a real measurement and was mistaken for one.
  //
  // So: a missing or non-executable binary is a hard error, and a non-zero exit
  // is reported on stderr rather than swallowed. Producing no bicliques remains
  // perfectly valid and is left to the caller to interpret.
  void runVnmextract(const std::string& graph_path, const std::string& out_prefix) const {
    const std::string exe = VNMEXTRACT_EXE;
    std::error_code ec;
    if (!std::filesystem::exists(exe, ec) || ec) {
      throw std::runtime_error(
          "BCCodec: vnmextract not found at '" + exe +
          "'. PDL-BC cannot be built. Ensure the dsextract external project is "
          "built (cmake --build <dir> --target dsextract).");
    }

    // <graph> format=1 shingle_size=1 <iters> <bcsizes> <out> <num_hashes>
    std::ostringstream cmd;
    cmd << exe << " " << graph_path << " 1 1 " << mining_.min_bicliques << " "
        << mining_.bcsizes << " " << out_prefix << " " << mining_.num_hashes
        << " >/dev/null 2>&1";
    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
      // Not fatal -- the caller falls back to verbatim sets -- but it must be
      // visible, because a silent failure looks exactly like "nothing to
      // compress" in the resulting index.
      std::cerr << "BCCodec: warning: vnmextract exited with status " << rc
                << " (graph=" << graph_path << "); continuing with whatever "
                   "biclique files it produced." << std::endl;
    }
  }

  // Parse `<prefix>-biclique-it-N.txt` files for N = 0, 1, 2, ... until
  // a missing file. Each line is `block_ids - doc_ids` (1-based block
  // node ids per drl's nodeIsBlock convention; 1-based doc ids).
  static void parseBicliqueFiles(
      const std::string& prefix,
      std::size_t n_doc,
      std::size_t n_slots,
      std::vector<std::vector<std::size_t>>& rule_docs,
      std::vector<std::vector<std::size_t>>& slot_rule_refs) {
    for (std::size_t it = 0;; ++it) {
      auto fn = prefix + "-biclique-it-" + std::to_string(it) + ".txt";
      std::ifstream in(fn);
      if (!in) break;
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        // Split on '-'. drl uses parseLine(*siter, before, '-', after).
        auto dash = line.find('-');
        if (dash == std::string::npos) continue;
        auto before = parseInts(line.substr(0, dash));
        auto after = parseInts(line.substr(dash + 1));
        if (before.empty() || after.empty()) continue;

        std::vector<std::size_t> docs;
        docs.reserve(after.size());
        for (auto v : after) {
          if (v == 0 || v > n_doc) continue;  // 1-based; sanity-clip.
          docs.push_back(v - 1);  // back to 0-based.
        }
        if (docs.empty()) continue;
        std::sort(docs.begin(), docs.end());
        docs.erase(std::unique(docs.begin(), docs.end()), docs.end());

        const std::size_t rule_idx = rule_docs.size();
        rule_docs.push_back(std::move(docs));

        for (auto v : before) {
          if (v <= n_doc) continue;  // not a block node.
          std::size_t slot = v - n_doc - 1;
          if (slot < n_slots) slot_rule_refs[slot].push_back(rule_idx);
        }
      }
    }
  }

  static std::vector<std::size_t> parseInts(const std::string& s) {
    std::vector<std::size_t> out;
    std::istringstream iss(s);
    std::size_t v;
    while (iss >> v) out.push_back(v);
    return out;
  }

  void packRules(const std::vector<std::vector<std::size_t>>& rule_docs, std::size_t n_doc) {
    std::size_t total = 0;
    for (const auto& r : rule_docs) total += r.size();
    sdsl::int_vector<> tmp_rules(total, 0);
    sdsl::bit_vector tmp_borders(std::max<std::size_t>(total, 1), 0);
    std::size_t pos = 0;
    for (const auto& r : rule_docs) {
      tmp_borders[pos] = 1;
      for (auto d : r) tmp_rules[pos++] = d;
    }
    sdsl::util::bit_compress(tmp_rules);
    rules_ = TIntVector(tmp_rules);
    rule_borders_ = TBitvector(tmp_borders);
    (void)n_doc;
  }

  void packBlocks(const std::vector<std::vector<std::size_t>>& slot_blocks) {
    std::size_t total = 0;
    for (const auto& b : slot_blocks) total += b.size();
    sdsl::int_vector<> tmp_blocks(std::max<std::size_t>(total, 1), 0);
    sdsl::bit_vector tmp_borders(std::max<std::size_t>(total, 1), 0);
    std::size_t pos = 0;
    for (const auto& b : slot_blocks) {
      tmp_borders[pos] = 1;
      for (auto v : b) tmp_blocks[pos++] = v;
    }
    if (total == 0) tmp_blocks.resize(0);
    sdsl::util::bit_compress(tmp_blocks);
    blocks_ = TIntVector(tmp_blocks);
    block_borders_ = TBitvector(tmp_borders);
  }

  // Not serialized: it only affects Build, never the built structure.
  MiningParams mining_;

  std::size_t n_doc_ = 0;
  std::size_t n_slots_ = 0;
  std::size_t n_rules_ = 0;
  TIntVector blocks_;
  TBitvector block_borders_;
  TSelect1 block_select_;
  TIntVector rules_;
  TBitvector rule_borders_;
  TSelect1 rule_select_;
};

static_assert(SetCodec<BCCodec<>>,
              "BCCodec<> must satisfy the SetCodec concept");

// collectSizes overloads — declared before any class-template that calls
// them unqualified (per the two-phase-lookup discipline established for
// grammar::CompactBPSLP etc. in include/dret/slp/compact_bp_slp.h). Each
// delegates to the codec's own GetSizeReport() and prepends the prefix.
template <typename TObjContainer, typename TPosContainer>
void collectSizes(SizeReport& out,
                  const PlainCodec<TObjContainer, TPosContainer>& codec,
                  const std::string& prefix = "") {
  for (const auto& f : codec.GetSizeReport()) append(out, prefix + f.name, f.bytes);
}

template <typename TSLP, typename TChunks>
void collectSizes(SizeReport& out,
                  const RPCodec<TSLP, TChunks>& codec,
                  const std::string& prefix = "") {
  for (const auto& f : codec.GetSizeReport()) append(out, prefix + f.name, f.bytes);
}

template <typename TBitvector, typename TIntVector>
void collectSizes(SizeReport& out,
                  const BCCodec<TBitvector, TIntVector>& codec,
                  const std::string& prefix = "") {
  for (const auto& f : codec.GetSizeReport()) append(out, prefix + f.name, f.bytes);
}

}  // namespace dret::pdl
