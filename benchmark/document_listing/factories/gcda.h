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
using SLP_Combined        = grammar::CombinedSLPWithUnitCover<>;

// Storage-parameterised typed-index template alias. The default TSLP matches
// the dret::gcda::DocListIdxGCDA default (SLP_Light = grammar::LightSLP<...>).
template <typename TStorage, typename TSLP = SLP_Light>
using Idx = dret::gcda::DocListIdxGCDA<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP>;

// Factory-path entry. Builds (constructs caches if missing, then loads)
// a typed index variant selected by `slp` and returns the type-erased
// DocListIndex<> handle plus its serialised byte size.
template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage,
     dret::Config& t_config,
     uint32_t t_block_size,
     float t_storing_factor,
     bench::axes::GCDASLPVariant t_slp) {
  std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t> result;

  auto build = [&](auto type_tag) {
    using TIndex = typename decltype(type_tag)::type;
    auto idx = std::make_shared<TIndex>(t_storage, t_block_size, t_storing_factor);
    construct(*idx, t_config);
    idx->load(t_config);
    result = {idx, sdsl::size_in_bytes(*idx)};
  };

  struct T_Light      { using type = Idx<TStorage, SLP_Light>; };
  struct T_CompactBP    { using type = Idx<TStorage, SLP_CompactBP>; };
  struct T_CompactLOUDS { using type = Idx<TStorage, SLP_CompactLOUDS>; };
  struct T_CSLP         { using type = Idx<TStorage, SLP_Combined>; };

  switch (t_slp) {
    case bench::axes::GCDASLPVariant::CompactBP:    build(T_CompactBP{});    break;
    case bench::axes::GCDASLPVariant::CompactLOUDS: build(T_CompactLOUDS{}); break;
    case bench::axes::GCDASLPVariant::Combined:         build(T_CSLP{});         break;
    case bench::axes::GCDASLPVariant::Light:
    default:                                        build(T_Light{});      break;
  }
  return result;
}

}  // namespace bench::factories::gcda
