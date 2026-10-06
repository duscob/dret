//
// GCDA factory header — type aliases + Make() for the
// dret::gcda::DocListIdxGCDA family.
//
// Storage-parametric type aliases let bm_query_doc_list (via the
// Factory<>'s shared ExternalGenericStorage) and bm_build_items (per-bench
// dret::GenericStorage) reuse one definition of the typed index. The Make()
// function dispatches the GCDASLPVariant axis to the right TSLP at compile
// time using the type-tag-lambda pattern, and runs construct() + load() so
// missing cache artefacts get built inline.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include <grammar/sampled_slp.h>
#include <grammar/slp.h>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_gcda_backend.h"
#include "dret/pdl/get_docs.h"
#include "dret/slp/combined_slp_with_unit_cover.h"
#include "dret/slp/compact_bp_slp.h"
#include "dret/slp/compact_louds_slp.h"

#include "../axes.h"

namespace bench::factories::gcda {

// TSLP type choices for the GCDA family. Storage-independent.
using SLP_Light     = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                          grammar::SampledSLP<>,
                                          grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
using SLP_CompactBP   = grammar::CompactBPSLP<>;
using SLP_CompactLOUDS = grammar::CompactLOUDSSLP<>;
// Bit-compressed (int_vector) rules, lengths, and leaves — matches light/compact.
// The all-default grammar::CombinedSLPWithUnitCover<> would store these as
// uncompressed std::vector<uint32_t> (~8x larger grammar base).
using SLP_Combined        = grammar::CombinedSLPWithUnitCover<
    grammar::CombinedSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
                         grammar::SampledSLP<>,
                         sdsl::int_vector<>>>;


// Document-set containers (GCDASetsCodec): Re-Pair compressed, the family
// default, or plain -- the sorted lists, bit-packed.
using Sets_RP = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                  true,
                                  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
using Sets_Plain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;

// Storage-parameterised typed-index template alias. The default TSLP matches
// the dret::gcda::DocListIdxGCDA default (SLP_Light = grammar::LightSLP<...>).
template <typename TStorage, typename TSLP = SLP_Light, typename TSets = Sets_RP>
using Idx = dret::gcda::DocListIdxGCDA<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP,
    TSets>;

// Factory-path entry. Builds (constructs caches if missing, then loads)
// a typed index variant selected by `slp` and returns the type-erased
// DocListIndex<> handle plus its serialised byte size.
template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage,
     dret::Config& t_config,
     uint32_t t_block_size,
     float t_storing_factor,
     bench::axes::GCDASLPVariant t_slp,
     bench::axes::GCDASetsCodec t_sets = bench::axes::GCDASetsCodec::RP) {
  std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t> result;

  auto build = [&](auto type_tag) {
    using TIndex = typename decltype(type_tag)::type;
    auto idx = std::make_shared<TIndex>(t_storage, t_block_size, t_storing_factor);
    construct(*idx, t_config);
    idx->load(t_config);
    result = {idx, sdsl::size_in_bytes(*idx)};
  };

  auto by_slp = [&]<typename TSets>() {
    struct T_Light        { using type = Idx<TStorage, SLP_Light, TSets>; };
    struct T_CompactBP    { using type = Idx<TStorage, SLP_CompactBP, TSets>; };
    struct T_CompactLOUDS { using type = Idx<TStorage, SLP_CompactLOUDS, TSets>; };
    struct T_CSLP         { using type = Idx<TStorage, SLP_Combined, TSets>; };
    switch (t_slp) {
      case bench::axes::GCDASLPVariant::CompactBP:    build(T_CompactBP{});    break;
      case bench::axes::GCDASLPVariant::CompactLOUDS: build(T_CompactLOUDS{}); break;
      case bench::axes::GCDASLPVariant::Combined:     build(T_CSLP{});         break;
      case bench::axes::GCDASLPVariant::Light:
      default:                                        build(T_Light{});        break;
    }
  };
  if (t_sets == bench::axes::GCDASetsCodec::Plain)
    by_slp.template operator()<Sets_Plain>();
  else
    by_slp.template operator()<Sets_RP>();
  return result;
}

// GCDA over a document-array backend: the sampled tree and node lists of the
// (block size, storing factor) cell, the ends of a range expanded by t_get_doc
// (DA, SA-Phi or RLCSA) instead of the grammar.
template <typename TStorage, typename TGetDocs, typename TSets>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
BuildBackend(TStorage t_storage, dret::Config& t_config, uint32_t t_block_size, float t_storing_factor) {
  using TIndex = dret::gcda::DocListIdxGCDABackend<TStorage, dret::Alphabet<>,
                                                   sri::RIndexCount<TStorage, dret::Alphabet<>>, TGetDocs, TSets>;
  auto idx = std::make_shared<TIndex>(t_storage, t_block_size, t_storing_factor);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <typename TStorage, typename TSets>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeBackendForSets(TStorage t_storage, dret::Config& t_config, uint32_t t_block_size, float t_storing_factor,
                   bench::axes::GetDocEnum t_get_doc) {
  constexpr auto kW = dret::Alphabet<>::int_width;
  switch (t_get_doc) {
    case bench::axes::GetDocEnum::SAPhiR:
      return BuildBackend<TStorage, dret::pdl::PDLGetDocsSAPhi_R<TStorage, kW>, TSets>(t_storage, t_config, t_block_size,
                                                                                       t_storing_factor);
    case bench::axes::GetDocEnum::RLCSA:
      return BuildBackend<TStorage, dret::pdl::PDLGetDocsRLCSA<TStorage, kW>, TSets>(t_storage, t_config, t_block_size,
                                                                                     t_storing_factor);
    case bench::axes::GetDocEnum::DA:
      return BuildBackend<TStorage, dret::pdl::PDLGetDocsDA<TStorage, kW>, TSets>(t_storage, t_config, t_block_size,
                                                                                  t_storing_factor);
    default:
      throw std::invalid_argument("GCDA backend: get_doc must be da, sa_phi_r or rlcsa");
  }
}

template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeBackend(TStorage t_storage,
            dret::Config& t_config,
            uint32_t t_block_size,
            float t_storing_factor,
            bench::axes::GetDocEnum t_get_doc,
            bench::axes::GCDASetsCodec t_sets = bench::axes::GCDASetsCodec::RP) {
  if (t_sets == bench::axes::GCDASetsCodec::Plain)
    return MakeBackendForSets<TStorage, Sets_Plain>(t_storage, t_config, t_block_size, t_storing_factor, t_get_doc);
  return MakeBackendForSets<TStorage, Sets_RP>(t_storage, t_config, t_block_size, t_storing_factor, t_get_doc);
}

}  // namespace bench::factories::gcda
