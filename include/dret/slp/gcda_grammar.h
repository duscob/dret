//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/21/26.
//
// Building the sampled grammars of GCDA: the sampled tree of a (block size, storing factor) cell over
// the CNF grammar of the DA, each representation of it (sampled-cached / LightSLP, sampled-ondemand,
// compact BP and LOUDS, differential), and the document lists of its sampled nodes in each codec. GCDA
// (doc_list/doc_list_gcda.h), the grammar backends (backend/grammar.h) and Brute-GCDA
// (slp/bare_grammar.h) build through these. The DA's grammar itself comes from slp/da_grammar.h.
//

#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/memory_management.hpp>

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "dret/cache_components.h"
#include "dret/config.h"
#include "dret/slp/combined_slp_with_unit_cover.h"
#include "dret/slp/compact_bp_slp.h"
#include "dret/slp/compact_louds_slp.h"
#include "dret/slp/da_grammar.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/slp/slp_components.h"

namespace dret::gcda {

// Key of the document lists of a cell's sampled nodes in a codec: bit-packed
// (plain) or grammar-compressed (Re-Pair). The raw lists every codec is built
// from are cached under the plain key in grammar::Chunks<>'s encoding.
template <typename TSets>
std::string NodeDocListsKey(const JSON& t_keys, const SampledTreeCell& t_cell) {
  constexpr bool kPlain = std::is_same_v<TSets, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>> ||
                          std::is_same_v<TSets, grammar::Chunks<>>;
  return CellKey(t_keys, kPlain ? conf::kDaNodeDocListsPlain : conf::kDaNodeDocListsRP, t_cell);
}

//~~~~~~~


template <bool kExpand, typename TChunks>
void construct(grammar::GCChunks<grammar::SLP<>, kExpand, TChunks>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename TSLP, bool kExpand>
void construct(grammar::GCChunks<TSLP, kExpand, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor);

// Plain (bit-packed, not grammar-compressed) document sets: the sorted lists
// themselves, one chunk per sampled node, as PDL's plain codec stores them.
inline void construct(grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>& t_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename TSLP, typename TSampledSLP, typename TLeavesContainer>
void construct(grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>& t_cslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename TSLP, typename TSampledSLP, typename TChunks>
void construct(grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

// Compact-grammar TSLP variants (from the grammar library; the dret-side
// shims at "dret/slp/compact_bp_slp.h" / "dret/slp/compact_louds_slp.h" re-export them as
// dret::CompactBPSLP / dret::CompactLOUDSSLP). These forward declarations are
// required so the unqualified `construct(slp, ...)` call inside the
// `construct(DocListIdxGCDA&, ...)` template below resolves at the
// template-definition point — ordinary lookup is fixed there, and ADL on
// `grammar::CompactBPSLP` searches `namespace grammar` only, never
// `dret::gcda` where the actual overload definitions live.
template <typename... Ts>
void construct(grammar::CompactBPSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename... Ts>
void construct(grammar::CompactLOUDSSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

// CombinedSLP-with-unit-cover variant. Same signature pattern as
// the compact-grammar variants above; same forward-decl placement reason
// (visibility for the dispatch template's unqualified `construct(slp, ...)`
// call below).
template <typename... Ts>
void construct(grammar::CombinedSLPWithUnitCover<Ts...>& t_wrapped,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

//~~~~~~~


// SLP-variant trait full definition. Lives below the construct(SLPType&, ...)
// forward declarations so its constructSLP body sees them via ordinary lookup
// at template-definition time. The primary template covers every
// GCDA-family SLP (grammar::LightSLP, CompactBPSLP, CompactLOUDSSLP,
// CombinedSLPWithUnitCover, ...): the 5-arg SLP-build path with a datafile
// parameter. The DifferentialLightSLP specialisation routes to the 4-arg
// DifferentialLightSLP::construct(). Both cache through CacheComponents.
template <typename TSLP>
struct GCDAVariantTraits {
  static void constructSLP(TSLP& slp,
                           Config& t_config,
                           const std::string& t_datafile,
                           uint32_t t_block_size,
                           float t_storing_factor) {
    construct(slp, t_config, t_datafile, t_block_size, t_storing_factor);
  }
};

template <typename... Ts>
struct GCDAVariantTraits<DifferentialLightSLP<Ts...>> {
  static void constructSLP(DifferentialLightSLP<Ts...>& slp,
                           Config& t_config,
                           const std::string& /*t_datafile_ignored*/,
                           uint32_t t_block_size,
                           float t_storing_factor) {
    // DifferentialLightSLP::construct (defined in
    // dret/slp/differential_light_slp.h) loads the DA from the config cache
    // itself, so no datafile parameter.
    construct(slp, t_config, t_block_size, t_storing_factor);
  }
};


//~~~~~~~


// The sampled tree of a (block size, storing factor) cell over the CNF grammar
// of the DA, with the document lists of its sampled nodes. Every representation
// derives from it; it is cached as components (the CNF grammar, the tree, its
// leaves), and the raw lists under the plain lists key, in grammar::Chunks<>'s
// encoding and bit-packed.
template <typename TSLP, typename TSampledSLP, typename TLeavesContainer>
void construct(grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>& t_cslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto components = CacheComponents(t_cslp, t_config.keys, cell);
  auto event = sdsl::memory_monitor::event(ComponentFile(components.back(), t_config));

  t_cslp = grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>(LoadOrBuildDaCnfGrammar(t_config, t_datafile));
  grammar::Chunks<> cslp_docs;

  grammar::AddSet add_set(cslp_docs);
  t_cslp.Compute(t_block_size,
                 add_set,
                 add_set,
                 grammar::MustBeSampled<decltype(cslp_docs)>(grammar::AreChildrenTooBig(cslp_docs, t_storing_factor)));

  StoreComponents(t_cslp, CacheComponents(t_cslp, t_config.keys, cell), t_config);

  const auto key_docs = CellKey(t_config.keys, kDaNodeDocListsPlain, cell);
  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);
}

// The raw sampled tree of a cell (see above): loaded from its components, or built.
inline grammar::CombinedSLP<> LoadOrBuildSampledTree(Config& t_config,
                                                     const std::string& t_datafile,
                                                     uint32_t t_block_size,
                                                     float t_storing_factor) {
  grammar::CombinedSLP<> cslp;
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  if (!LoadComponents(cslp, CacheComponents(cslp, t_config.keys, cell), t_config))
    construct(cslp, t_config, t_datafile, t_block_size, t_storing_factor);
  return cslp;
}

// Sampled-cached: the DA grammar (not in CNF) with, for each leaf of the cell's
// sampled tree, its cover in that grammar.
template <typename TSLP, typename TSampledSLP, typename TChunks>
void construct(grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto cslp = LoadOrBuildSampledTree(t_config, t_datafile, t_block_size, t_storing_factor);

  // Collection-level, not cell-level: parsed once and reused by every
  // (bs,sf) cell. See LoadOrBuildDaSlp above.
  grammar::SLP<> slp;
  std::vector<std::size_t> compact_seq;
  LoadOrBuildDaSlp(t_config, t_datafile, slp, compact_seq);

  auto event = sdsl::memory_monitor::event(ComponentFile(CacheComponents(t_lslp, t_config.keys, cell).back(), t_config));
  grammar::LightSLP<> lslp;
  lslp.Compute(slp, compact_seq, cslp);

  auto bit_compress = [](auto& _v) {
    sdsl::util::bit_compress(_v);
  };
  t_lslp = grammar::LightSLP<TSLP, TSampledSLP, TChunks>(lslp, bit_compress, bit_compress, bit_compress, bit_compress);
  StoreComponents(t_lslp, CacheComponents(t_lslp, t_config.keys, cell), t_config);
}

//~~~~~~~


// Shared helper for the compact-grammar construct overloads. Both BP and
// LOUDS variants derive the compact class from the cell's raw sampled tree via
// Compute(): a compact encoding of the CNF grammar (one per collection), the
// sampled leaves numbered in it, and the shared sampled tree.
template <typename TCompactSLP>
void constructCompactCommon(TCompactSLP& t_compact,
                            Config& t_config,
                            const std::string& t_datafile,
                            uint32_t t_block_size,
                            float t_storing_factor) {
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto cslp = LoadOrBuildSampledTree(t_config, t_datafile, t_block_size, t_storing_factor);

  auto event = sdsl::memory_monitor::event(ComponentFile(CacheComponents(t_compact, t_config.keys, cell)[1], t_config));
  t_compact.Compute(cslp);
  StoreComponents(t_compact, CacheComponents(t_compact, t_config.keys, cell), t_config);
}

template <typename... Ts>
void construct(grammar::CompactBPSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  constructCompactCommon(t_compact, t_config, t_datafile, t_block_size, t_storing_factor);
}

template <typename... Ts>
void construct(grammar::CompactLOUDSSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  constructCompactCommon(t_compact, t_config, t_datafile, t_block_size, t_storing_factor);
}

// Sampled-ondemand: the cell's raw sampled tree converted into the wrapper's
// (possibly bit-compressed) CombinedSLP base -- applying `bit_compress` to the
// SLP rules, lengths, and leaves.
template <typename... Ts>
void construct(grammar::CombinedSLPWithUnitCover<Ts...>& t_wrapped,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using BaseCSLP = typename grammar::CombinedSLPWithUnitCover<Ts...>::Base;
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto raw = LoadOrBuildSampledTree(t_config, t_datafile, t_block_size, t_storing_factor);

  auto event = sdsl::memory_monitor::event(ComponentFile(CacheComponents(t_wrapped, t_config.keys, cell).back(), t_config));

  // The production SLP_Combined uses sdsl::int_vector<> containers; the action is a
  // no-op for any non-int_vector container (e.g. the default std::vector base).
  auto bit_compress = [](auto& v) {
    if constexpr (std::is_same_v<std::decay_t<decltype(v)>, sdsl::int_vector<>>)
      sdsl::util::bit_compress(v);
  };
  BaseCSLP cslp(raw, bit_compress, bit_compress, bit_compress);

  static_cast<BaseCSLP&>(t_wrapped) = cslp;
  StoreComponents(t_wrapped, CacheComponents(t_wrapped, t_config.keys, cell), t_config);
}

//~~~~~~~


// The raw lists of a cell, cached while its sampled tree was built.
inline grammar::Chunks<> LoadRawNodeDocLists(Config& t_config, uint32_t t_block_size, float t_storing_factor) {
  const auto key_docs = NodeDocListsKey<grammar::Chunks<>>(t_config.keys, {t_block_size, t_storing_factor});
  grammar::Chunks<> sets;
  if (!sdsl::load_from_cache(sets, key_docs, t_config, true))
    throw std::runtime_error("GCDA document lists: no cached lists under " + key_docs +
                             "; build the sampled tree for this (block size, storing factor) first");
  return sets;
}

template <bool kExpand, typename TChunks>
void construct(grammar::GCChunks<grammar::SLP<>, kExpand, TChunks>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor) {
  const auto slp_sets = LoadRawNodeDocLists(t_config, t_block_size, t_storing_factor);
  const auto& objs = slp_sets.GetObjects();
  grammar::RePairEncoder<false> encoder_nslp;
  t_slp_sets.Compute(objs.begin(), objs.end(), slp_sets, encoder_nslp);
  sdsl::store_to_cache(t_slp_sets, NodeDocListsKey<std::decay_t<decltype(t_slp_sets)>>(t_config.keys,
                                                                                        {t_block_size, t_storing_factor}),
                       t_config, true);
}

//~~~~~~~


template <typename TSLP, bool kExpand>
void construct(grammar::GCChunks<TSLP, kExpand, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor) {
  const auto key_docs = NodeDocListsKey<std::decay_t<decltype(t_slp_sets)>>(t_config.keys,
                                                                            {t_block_size, t_storing_factor});

  grammar::GCChunks<grammar::SLP<>> slp_sets;
  if (!sdsl::cache_file_exists<decltype(slp_sets)>(key_docs, t_config)) {
    construct(slp_sets, t_config, t_block_size, t_storing_factor);
  } else {
    sdsl::load_from_cache(slp_sets, key_docs, t_config, true);
  }

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };

  t_slp_sets =
      std::remove_reference_t<decltype(t_slp_sets)>(slp_sets, bit_compress, bit_compress, bit_compress, bit_compress);
  sdsl::store_to_cache(t_slp_sets, key_docs, t_config, true);
}

//~~~~~~~


// Plain (bit-packed, not grammar-compressed) document lists: the sorted lists
// themselves, one chunk per sampled node, as PDL's plain codec stores them.
// Cached as a side effect of building the sampled tree; GCDA's construct
// reaches this only when the bit-packed copy is missing.
inline void construct(grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>& t_sets,
                      Config& t_config,
                      uint32_t t_block_size,
                      float t_storing_factor) {
  const auto sets = LoadRawNodeDocLists(t_config, t_block_size, t_storing_factor);
  auto bit_compress = [](sdsl::int_vector<>& v) { sdsl::util::bit_compress(v); };
  t_sets = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>(sets, bit_compress, bit_compress);
  sdsl::store_to_cache(t_sets, NodeDocListsKey<grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>(
                                   t_config.keys, {t_block_size, t_storing_factor}),
                       t_config, true);
}

//~~~~~~~


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

}  // namespace dret::gcda
