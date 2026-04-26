//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// dret-side shim for grammar::CompactLOUDSSLP. The class itself lives
// upstream (added in Phase A.0.3); this header only:
//   - re-exports it as `dret::CompactLOUDSSLP` for downstream wiring,
//   - declares the dret-side `collectSizes` overload (Task 7),
//   - declares a dedicated `dret::ExpandSLP` overload (Task 6) that
//     bypasses the default recursive-descent path (which would walk
//     LOUDS' O(left subtree) `operator[]` for every internal node) and
//     instead delegates to `grammar::ExpandCompactSLPForward`, which
//     uses excess tracking for an O(N) total scan of the requested range.
//
// Must be `#include`-d at the top of `doc_list_sampled_tree_gcda.h` so
// `DocListIdxGCDA::getDocs` (unqualified `ExpandSLP(*slp_, ...)`) and
// `GetSizeReport` (unqualified `collectSizes(*slp_, ...)`) resolve at
// the template-definition point.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>

#include <grammar/compact_louds_slp.h>
#include <grammar/slp_compact.h>
#include <sdsl/util.hpp>

#include "size_report.h"

namespace dret {

template<typename... Ts>
using CompactLOUDSSLP = grammar::CompactLOUDSSLP<Ts...>;

// ---- size report ----

template<typename... Ts>
void collectSizes(SizeReport& out,
                  const grammar::CompactLOUDSSLP<Ts...>& slp,
                  const std::string& prefix = "") {
  append(out, prefix + "sigma",          sizeof(std::size_t));
  append(out, prefix + "bv_tree",        sdsl::size_in_bytes(slp.BvTree()));
  append(out, prefix + "leaf_rank",      sdsl::size_in_bytes(slp.LeafRank()));
  append(out, prefix + "compact_leaves", sdsl::size_in_bytes(slp.CompactLeaves()));
  append(out, prefix + "sampled_leaves", sdsl::size_in_bytes(slp.SampledLeaves()));
  append(out, prefix + "sampled_slp",
         sdsl::size_in_bytes(static_cast<const grammar::SampledSLP<>&>(slp)));
}

// ---- ExpandSLP overload ----
//
// Dedicated forward-DFS path for LOUDS. The default `dret::ExpandSLP`
// (slp_tools.h:132) decomposes the [bp, ep) range into per-leaf calls of
// `ExpandSLPFromFront`/`ExpandSLPFromBack`, which call `slp.Cover(leaf)`
// then recurse via `slp[var]`. For LOUDS each `slp[var]` is O(left
// subtree), so the recursive path is O(N log n) for N terminals on a
// balanced grammar. `grammar::ExpandCompactSLPForward` walks the bit
// vector linearly with excess tracking — O(N) total — and stops when its
// internal `t_length` counter reaches zero (no exception needed).
//
// Algorithm: iterate through the sampled-tree leaves intersecting [bp, ep).
// For each leaf, expand its `Map` variable forward by either the full leaf
// length or the remaining `ep - bp`, whichever is smaller. The first leaf
// may have a non-zero `skip_front` (when bp lies mid-leaf); in that case
// we wrap the user's `report` callback to drop the first `skip_front`
// invocations before forwarding the next ones.
//
// `grammar::ExpandCompactSLPForward` takes a non-const lvalue `TReport &`,
// so the wrapper must be a named local lambda (named-lambda pattern).
template<typename... Ts, typename Report>
void ExpandSLP(const grammar::CompactLOUDSSLP<Ts...>& slp,
               std::size_t bp,
               std::size_t ep,
               Report& report) {
  if (bp >= ep) return;

  auto leaf = slp.Leaf(bp);
  auto pos = slp.Position(leaf);
  std::size_t skip = bp - pos;  // 0 if front-aligned, else mid-leaf offset

  while (bp < ep) {
    auto next_pos = slp.Position(leaf + 1);
    auto var = slp.Map(leaf);
    auto target = std::min(ep, next_pos);
    auto take = target - bp;            // terminals to emit from this leaf
    auto length = skip + take;          // expansion length including skip prefix

    if (skip == 0) {
      grammar::ExpandCompactSLPForward(slp.BvTree(),
                                        slp.Sigma(),
                                        slp.CompactLeaves(),
                                        slp.LeafRank(),
                                        var,
                                        length,
                                        report);
    } else {
      std::size_t remaining_skip = skip;
      auto wrapped = [&report, &remaining_skip](const auto& v) {
        if (remaining_skip > 0) { --remaining_skip; return; }
        report(v);
      };
      grammar::ExpandCompactSLPForward(slp.BvTree(),
                                        slp.Sigma(),
                                        slp.CompactLeaves(),
                                        slp.LeafRank(),
                                        var,
                                        length,
                                        wrapped);
    }

    bp = target;
    ++leaf;
    pos = next_pos;
    skip = 0;  // only the first leaf has a skip prefix
  }
}

}  // namespace dret
