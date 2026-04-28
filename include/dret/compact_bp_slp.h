//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// dret-side shim for grammar::CompactBPSLP. The class itself lives upstream
// (added in Phase A.0.3); this header only:
//   - re-exports it as `dret::CompactBPSLP` for downstream wiring,
//   - declares the dret-side `collectSizes` overload so two-phase lookup at
//     `DocListIdxGCDA::GetSizeReport`'s definition point can find it.
//
// Must be `#include`-d at the top of `doc_list_sampled_tree_gcda.h` (before
// the `DocListIdxGCDA` class body) so the size-report overload is visible to
// ordinary lookup at the template-definition point. ADL on
// `grammar::CompactBPSLP` searches `namespace grammar` only — never `dret` —
// so without this declaration `collectSizes(*slp_, ...)` would fail.
//

#pragma once

#include <string>

#include <grammar/compact_bp_slp.h>
#include <sdsl/util.hpp>

#include "size_report.h"

namespace dret {

template<typename... Ts>
using CompactBPSLP = grammar::CompactBPSLP<Ts...>;

template<typename... Ts>
void collectSizes(SizeReport& out,
                  const grammar::CompactBPSLP<Ts...>& slp,
                  const std::string& prefix = "") {
  append(out, prefix + "sigma",          sizeof(std::size_t));
  append(out, prefix + "bv_tree",        sdsl::size_in_bytes(slp.BvTree()));
  append(out, prefix + "bp_support",     sdsl::size_in_bytes(slp.BpSupport()));
  append(out, prefix + "bv_leaf_marks",  sdsl::size_in_bytes(slp.BvLeafMarks()));
  append(out, prefix + "bv_leaf_marks_rank",
         sdsl::size_in_bytes(slp.BvLeafMarksRank()));
  append(out, prefix + "compact_leaves", sdsl::size_in_bytes(slp.CompactLeaves()));
  append(out, prefix + "sampled_leaves", sdsl::size_in_bytes(slp.SampledLeaves()));
  append(out, prefix + "sampled_slp",
         sdsl::size_in_bytes(static_cast<const grammar::SampledSLP<>&>(slp)));
}

}  // namespace dret
