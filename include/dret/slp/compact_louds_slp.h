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
// Must be `#include`-d at the top of `doc_list/doc_list_gcda.h` so
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

#include "dret/size_report.h"

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
// may have a non-zero `skip_front` (when bp lies mid-leaf). It is expanded
// from its nearer end -- backward when the range runs to the leaf's end or
// lies closer to it -- and the user's `report` callback is wrapped to drop
// the skipped invocations before forwarding the next ones. Backward reports
// the leaf's part in reverse order; ExpandSLPUntil keeps the forward order.
//
// `grammar::ExpandCompactSLPForward` takes a non-const lvalue `TReport &`,
// so the wrapper must be a named local lambda (named-lambda pattern).
// kNearerEnd = false keeps every leaf forward, in position order (ExpandSLPUntil).
template<bool kNearerEnd, typename... Ts, typename Report>
void ExpandLOUDSLeaves(const grammar::CompactLOUDSSLP<Ts...>& slp,
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
      // Reach bp from the nearer end of the leaf: LOUDS stores no span lengths, so the skipped
      // positions are expanded and dropped, and from the far end they would cost up to the leaf.
      std::size_t skip_back = next_pos - target;  // 0 when the range runs to the end of the leaf
      bool backward = kNearerEnd && skip_back < skip;
      std::size_t remaining_skip = backward ? skip_back : skip;
      auto wrapped = [&report, &remaining_skip](const auto& v) {
        if (remaining_skip > 0) { --remaining_skip; return; }
        report(v);
      };
      if (backward) {
        grammar::ExpandCompactSLPBackward(slp.BvTree(),
                                           slp.Sigma(),
                                           slp.CompactLeaves(),
                                           slp.LeafRank(),
                                           var,
                                           skip_back + take,
                                           wrapped);
      } else {
        grammar::ExpandCompactSLPForward(slp.BvTree(),
                                          slp.Sigma(),
                                          slp.CompactLeaves(),
                                          slp.LeafRank(),
                                          var,
                                          length,
                                          wrapped);
      }
    }

    bp = target;
    ++leaf;
    pos = next_pos;
    skip = 0;  // only the first leaf has a skip prefix
  }
}

template<typename... Ts, typename Report>
void ExpandSLP(const grammar::CompactLOUDSSLP<Ts...>& slp,
               std::size_t bp,
               std::size_t ep,
               Report& report) {
  ExpandLOUDSLeaves<true>(slp, bp, ep, report);
}

// ---- ExpandSLPUntil overload ----
//
// The forward-only expansion is in position order, but
// grammar::ExpandCompactSLPForward offers no way to stop, so this keeps the
// order and drops everything after the report asks to stop. Correct, not
// faster: LOUDS still expands the whole range.
template<typename... Ts, typename Report>
void ExpandSLPUntil(const grammar::CompactLOUDSSLP<Ts...>& slp,
                    std::size_t bp,
                    std::size_t ep,
                    Report& report) {
  bool going = true;
  auto filtered = [&report, &going](auto v) {
    if (going)
      going = report(v);
  };
  ExpandLOUDSLeaves<false>(slp, bp, ep, filtered);
}

}  // namespace dret
