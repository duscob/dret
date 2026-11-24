//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/4/26.
//
// Cache components of the grammars over the document array (DA): which range
// of each grammar class's serialization is which component, and under which
// key. Every index that contains a component reads the same file: the CNF
// grammar is shared by GCDA's on-demand representation, GCDA-nolists and the
// sampled-tree builds; the sampled tree of a (block size, storing factor) cell
// by every GCDA representation and GCDA-differential; and so on.
//
// The grammar classes serialize their bases first, in order (grammar, sampled
// tree, then the representation's own data; the compact BP / LOUDS classes put
// their grammar encoding first and the sampled tree last), so each component is
// one contiguous range and the classes' formats are unchanged.
//

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include <grammar/combined_slp_with_unit_cover.h>
#include <grammar/compact_bp_slp.h>
#include <grammar/compact_louds_slp.h>
#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>

#include "dret/cache_components.h"
#include "dret/config.h"
#include "dret/construct_base.h"

namespace dret {

// The parameters of a sampled tree: its block size and storing factor.
struct SampledTreeCell {
  uint32_t block_size = 512;
  float storing_factor = 4;
};

inline std::string KeyName(const JSON& t_keys, std::string_view t_name) {
  return t_keys[t_name].get<std::string>();
}

// Key of a component that depends on the sampled-tree cell.
inline std::string CellKey(const JSON& t_keys, std::string_view t_name, const SampledTreeCell& t_cell) {
  return PrefixedKey(t_keys, conf::kBlkSf, KeyName(t_keys, t_name), t_cell.block_size,
                     StoringFactorTerm(t_cell.storing_factor));
}

//~~~~~~~  The CNF grammar alone (GCDA-nolists, plain)  ~~~~~~~

template <typename TVars, typename TLengths>
std::vector<Component> CacheComponents(const grammar::SLP<TVars, TLengths>& t_slp, const JSON& t_keys) {
  return {{KeyName(t_keys, conf::kDaCnfGrammar), TypeHash<grammar::SLP<TVars, TLengths>>(), SerializedSize(t_slp)}};
}

//~~~~~~~  CNF grammar + sampled tree + leaves (sampled-ondemand)  ~~~~~~~

template <typename TSLP, typename TSampledSLP, typename TLeaves>
std::vector<Component> CacheComponents(const grammar::CombinedSLP<TSLP, TSampledSLP, TLeaves>& t_cslp,
                                       const JSON& t_keys,
                                       const SampledTreeCell& t_cell) {
  const auto grammar = SerializedSize(static_cast<const TSLP&>(t_cslp));
  const auto tree = SerializedSize(static_cast<const TSampledSLP&>(t_cslp));
  return {{KeyName(t_keys, conf::kDaCnfGrammar), TypeHash<TSLP>(), grammar},
          {CellKey(t_keys, conf::kDaSampledTree, t_cell), TypeHash<TSampledSLP>(), tree},
          {CellKey(t_keys, conf::kDaSampledLeaves, t_cell), TypeHash<TLeaves>(),
           SerializedSize(t_cslp) - grammar - tree}};
}

// The unit-cover wrapper adds no data: it is stored as its CombinedSLP.
template <typename TCombinedSLP>
std::vector<Component> CacheComponents(const grammar::CombinedSLPWithUnitCover<TCombinedSLP>& t_wrapped,
                                       const JSON& t_keys,
                                       const SampledTreeCell& t_cell) {
  return CacheComponents(static_cast<const TCombinedSLP&>(t_wrapped), t_keys, t_cell);
}

//~~~~~~~  DA grammar + sampled tree + leaf covers (sampled-cached)  ~~~~~~~

template <typename TSLP, typename TSampledSLP, typename TChunks>
std::vector<Component> CacheComponents(const grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
                                       const JSON& t_keys,
                                       const SampledTreeCell& t_cell) {
  const auto grammar = SerializedSize(static_cast<const TSLP&>(t_lslp));
  const auto tree = SerializedSize(static_cast<const TSampledSLP&>(t_lslp));
  return {{KeyName(t_keys, conf::kDaGrammar), TypeHash<TSLP>(), grammar},
          {CellKey(t_keys, conf::kDaSampledTree, t_cell), TypeHash<TSampledSLP>(), tree},
          {CellKey(t_keys, conf::kDaLeafCovers, t_cell), TypeHash<TChunks>(), SerializedSize(t_lslp) - grammar - tree}};
}

//~~~~~~~  Compact BP / LOUDS: grammar encoding + leaves in it + sampled tree  ~~~~~~~

namespace internal {

template <typename TCompact>
std::vector<Component> CompactComponents(const TCompact& t_compact,
                                         const JSON& t_keys,
                                         const SampledTreeCell& t_cell,
                                         std::string_view t_grammar,
                                         std::string_view t_leaves) {
  using TSampled = typename TCompact::SampledBase;
  const auto leaves = SerializedSize(t_compact.SampledLeaves());
  const auto tree = SerializedSize(static_cast<const TSampled&>(t_compact));
  // The grammar encoding and the leaves are numbered in it, so both carry the
  // compact class's type.
  return {{KeyName(t_keys, t_grammar), TypeHash<TCompact>(), SerializedSize(t_compact) - leaves - tree},
          {CellKey(t_keys, t_leaves, t_cell), TypeHash<TCompact>(), leaves},
          {CellKey(t_keys, conf::kDaSampledTree, t_cell), TypeHash<TSampled>(), tree}};
}

}  // namespace internal

template <typename... Ts>
std::vector<Component> CacheComponents(const grammar::CompactBPSLP<Ts...>& t_compact,
                                       const JSON& t_keys,
                                       const SampledTreeCell& t_cell) {
  return internal::CompactComponents(t_compact, t_keys, t_cell, conf::kDaCnfGrammarBP, conf::kDaSampledLeavesBP);
}

template <typename... Ts>
std::vector<Component> CacheComponents(const grammar::CompactLOUDSSLP<Ts...>& t_compact,
                                       const JSON& t_keys,
                                       const SampledTreeCell& t_cell) {
  return internal::CompactComponents(t_compact, t_keys, t_cell, conf::kDaCnfGrammarLOUDS,
                                     conf::kDaSampledLeavesLOUDS);
}

}  // namespace dret
