//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/11/26.
//
// PDL raw-range get-doc policy: the thin adapter PDL indexes use to read
// document ids for SA positions that the sparse suffix tree did not
// cover. The interface deliberately matches the existing
// dret::rmq::GetDoc* family (rmq_get_doc_policies.h) so PDLRawRangePolicy
// can wrap any of them without duplicating the SA-range iteration code.
// Tasks 20-22 wire the three concrete backings (DA, GCDA-SLP, DGCDA) by
// instantiating PDLRawRangePolicy<rmq::GetDocDA<>> /
// PDLRawRangePolicy<rmq::GetDocSLP<>> / PDLRawRangePolicy<rmq::GetDocDSLP<>>.
//

#pragma once

#include <cstddef>
#include <utility>

#include "../config.h"
#include "../rmq_get_doc_policies.h"

namespace dret::pdl {

// Wraps a raw-range source (DA, grammar-compressed DA, differential
// grammar-compressed DA — anything exposing operator()(i) and
// operator()(b, e, report)) and re-exposes it under the names
// PDL queries use. The wrapped type still owns its storage / cache
// loading; the adapter doesn't store any DA-like state of its own.
template <typename TInner>
class PDLRawRangePolicy {
 public:
  using Inner = TInner;

  PDLRawRangePolicy() = default;
  explicit PDLRawRangePolicy(TInner t_inner) : inner_(std::move(t_inner)) {}

  // Half-open SA-range expansion: invokes t_report(doc) for every
  // position in [t_sp, t_ep). Doc-id remapping (1-based to 0-based,
  // sentinel handling, etc.) lives inside TInner — the adapter just
  // forwards.
  template <typename TReport>
  void getDocs(std::size_t t_sp, std::size_t t_ep, TReport& t_report) const {
    inner_(t_sp, t_ep, t_report);
  }

  // Single-position lookup. Convenience for callers that already know
  // they want one position and would otherwise wrap a one-element loop.
  std::size_t getDoc(std::size_t t_i) const {
    return inner_(t_i);
  }

  TInner& inner() { return inner_; }
  const TInner& inner() const { return inner_; }

 private:
  TInner inner_{};
};

// Concrete instantiations wiring the three rmq::GetDoc* backings into
// PDLRawRangePolicy. Defaults mirror the underlying types so the typed
// cache entries are shared with the existing RMQ document-listing
// indexes (any divergence would cost us a duplicate cache file).
//
// Tasks 20 / 21 / 22 of docs/pdl_indexes_tasks.md.

// Task 20 — DA-backed.
template <typename TStorage = GenericStorage, uint8_t t_width = 8>
using PDLGetDocsDA = PDLRawRangePolicy<rmq::GetDocDA<TStorage, t_width>>;

// Task 21 — GCDA / LightSLP-backed.
template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>>
using PDLGetDocsSLP = PDLRawRangePolicy<rmq::GetDocSLP<TStorage, t_width, TSLP>>;

// Task 22 — DGCDA / DifferentialLightSLP-backed.
template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TDSLP = DifferentialLightSLP<>>
using PDLGetDocsDSLP = PDLRawRangePolicy<rmq::GetDocDSLP<TStorage, t_width, TDSLP>>;

// Free-function construct() that forwards to the inner rmq::GetDoc*
// construct() so PDL indexes can build their raw get-doc cache through
// the project's normal construction API. DA-backed is a no-op (the DA
// is already built upstream); SLP / DSLP-backed builds the
// corresponding LightSLP / DifferentialLightSLP cache entries.
template <typename TInner>
void construct(PDLRawRangePolicy<TInner>& t_policy, Config& t_config) {
  construct(t_policy.inner(), t_config);
}

}  // namespace dret::pdl
