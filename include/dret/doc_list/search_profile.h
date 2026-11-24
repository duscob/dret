//
// Optional, zero-overhead-by-default phase profiler for the sampled-tree
// document-listing Search() path (DLSampledTreeScheme).
//
// Enabled only when the translation unit is compiled with
// -DDRET_DOC_LIST_PROFILE; otherwise this header defines nothing and the
// instrumentation in doc_list_sampled_tree_base.h compiles out entirely.
//
// Purpose: split a query's time into its four phases so the combine-vs-expand
// trade-off behind the block-size U-shape (Page vs Revision) can be measured
// directly, together with per-query structural counts (covers, boundary
// positions, distinct documents). See the journal experiments discussion.
//

#pragma once

#ifdef DRET_DOC_LIST_PROFILE

#include <cstdint>

namespace dret {

// Single-threaded accumulator (the benchmark runs one thread). Times are in
// nanoseconds, summed over every Search() call since the last reset().
struct SearchProfile {
  std::uint64_t ns_count = 0;    // count(): map the pattern to the SA range
  std::uint64_t ns_cover = 0;    // computeCoverFull(): locate the in-range covers
  std::uint64_t ns_expand = 0;   // getDocs() over the boundary + its sort/unique
  std::uint64_t ns_combine = 0;  // merge_sets_(): union the precomputed cover sets

  std::uint64_t n_queries = 0;       // Search() calls
  std::uint64_t n_nodes = 0;         // covers combined (summed over queries)
  std::uint64_t n_raw_positions = 0; // boundary positions expanded (summed)
  std::uint64_t n_docs = 0;          // distinct documents returned (summed)

  void reset() { *this = SearchProfile{}; }
};

// Process-wide instance; reset by the benchmark before each timed sweep point.
inline SearchProfile& search_profile() {
  static SearchProfile p;
  return p;
}

}  // namespace dret

#endif  // DRET_DOC_LIST_PROFILE
