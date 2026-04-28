//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// dret-side shim for grammar::CombinedSLPWithUnitCover. The class itself
// lives upstream; this header re-exports it as
// `dret::CombinedSLPWithUnitCover` for downstream wiring and declares the
// dret-side `collectSizes` overload so two-phase lookup at
// `DocListIdxGCDA::GetSizeReport`'s definition point can find it.
//
// Must be `#include`-d at the top of `doc_list_sampled_tree_gcda.h`
// (before the `DocListIdxGCDA` class body).
//

#pragma once

#include <string>

#include <grammar/combined_slp_with_unit_cover.h>
#include <grammar/sampled_slp.h>
#include <sdsl/util.hpp>

#include "size_report.h"

namespace dret {

template<typename... Ts>
using CombinedSLPWithUnitCover = grammar::CombinedSLPWithUnitCover<Ts...>;

template<typename... Ts>
void collectSizes(SizeReport& out,
                  const grammar::CombinedSLPWithUnitCover<Ts...>& slp,
                  const std::string& prefix = "") {
  // The wrapper has no data of its own — all storage is in the inherited
  // grammar::CombinedSLP base, which is itself a SLP base + a SampledSLP
  // base + a leaves vector. Report them as the underlying CSLP would, so
  // the GCDA-CSLP size report stays comparable to GCDA-Default's CSLP
  // scaffold.
  using BaseCSLP = typename grammar::CombinedSLPWithUnitCover<Ts...>::Base;
  const BaseCSLP& base = static_cast<const BaseCSLP&>(slp);
  append(out, prefix + "base_slp", sdsl::size_in_bytes(static_cast<const grammar::SLP<>&>(base)));
  append(out, prefix + "leaves",   sdsl::size_in_bytes(base.GetLeaves()));
  append(out, prefix + "sampled_slp",
         sdsl::size_in_bytes(static_cast<const grammar::SampledSLP<>&>(base)));
}

}  // namespace dret
