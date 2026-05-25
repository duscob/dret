//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// dret-side shim for grammar::CombinedSLPWithUnitCover. The class itself
// lives upstream; this header re-exports it as
// `dret::CombinedSLPWithUnitCover` for downstream wiring and declares the
// dret-side `collectSizes` overload so two-phase lookup at
// `DocListIdxGCDA::GetSizeReport`'s definition point can find it.
//
// Must be `#include`-d at the top of `doc_list/doc_list_gcda.h`
// (before the `DocListIdxGCDA` class body).
//

#pragma once

#include <string>

#include <grammar/combined_slp_with_unit_cover.h>
#include <grammar/sampled_slp.h>
#include <sdsl/util.hpp>

#include "dret/size_report.h"

namespace dret {

template<typename... Ts>
using CombinedSLPWithUnitCover = grammar::CombinedSLPWithUnitCover<Ts...>;

// Deduces the actual SLP / SampledSLP / leaves container types from the wrapped
// grammar::CombinedSLP so the report reflects whatever containers the variant
// uses (e.g. bit-compressed sdsl::int_vector<>), not a hardcoded grammar::SLP<>.
template<typename TSLP, typename TSampledSLP, typename TLeaves>
void collectSizes(SizeReport& out,
                  const grammar::CombinedSLPWithUnitCover<
                      grammar::CombinedSLP<TSLP, TSampledSLP, TLeaves>>& slp,
                  const std::string& prefix = "") {
  // The wrapper has no data of its own — all storage is in the inherited
  // grammar::CombinedSLP base: a SLP base + a SampledSLP base + a leaves vector.
  using BaseCSLP = grammar::CombinedSLP<TSLP, TSampledSLP, TLeaves>;
  const BaseCSLP& base = static_cast<const BaseCSLP&>(slp);
  append(out, prefix + "base_slp",   sdsl::size_in_bytes(static_cast<const TSLP&>(base)));
  append(out, prefix + "leaves",     sdsl::size_in_bytes(base.GetLeaves()));
  append(out, prefix + "sampled_slp", sdsl::size_in_bytes(static_cast<const TSampledSLP&>(base)));
}

}  // namespace dret
