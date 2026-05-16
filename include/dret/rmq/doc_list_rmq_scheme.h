//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Stateless LEFTMOST-RMQ traversal and supporting hook classes for the
// document-listing RMQ indexes (SADA / ILCP / CILCP). Kept free of grammar /
// SLP dependencies so it can be included cheaply in any listing core.
//
// This header is intentionally *not* used by doc_freq/doc_freq_rmq.h — the
// frequency-index path continues to use GetExtremeOccurrencesRMQ from
// doc_index_rmq_common.h.
//

#pragma once

#include <cstddef>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

namespace dret::rmq {

// Recursive LEFTMOST RMQ traversal over half-open [bp, ep).
//
// Callbacks:
//   get_doc(k)     -> std::size_t   — doc id at SA position k
//   is_reported(k, d) -> bool       — true if doc d already reported
//   report(k, d)   -> void          — mark d as reported and emit it
//
// For SADA: pass MarkedReported as is_reported; a lambda that calls
//   mr.mark(d) + user_report(d) as report.
// For ILCP: report also fans out the rest of the run; see IlcpLikeCore.
template <typename TRMQ, typename TGetDoc, typename TIsReported, typename TReport>
void ListDocsRMQScheme(std::size_t bp,
                       std::size_t ep,
                       const TRMQ& rmq,
                       const TGetDoc& get_doc,
                       const TIsReported& is_reported,
                       TReport& report) {
  if (bp >= ep)
    return;
  const auto k = rmq(bp, ep - 1);
  const auto d = get_doc(k);
  if (!is_reported(k, d)) {
    report(k, d);
    ListDocsRMQScheme(bp, k, rmq, get_doc, is_reported, report);
    ListDocsRMQScheme(k + 1, ep, rmq, get_doc, is_reported, report);
  }
}

//~~~~~~~


// Tracks which documents have been reported in the current query.
// Serves as the is_reported predicate; marking is done explicitly via mark().
template <typename TBV = sdsl::bit_vector>
class MarkedReported {
 public:
  explicit MarkedReported(std::size_t n_doc) : marked_(n_doc, 0) {}

  // is_reported predicate: doc d >= n_doc counts as "already reported" so
  // sentinel suffixes (from the terminator document) are silently skipped.
  bool operator()(std::size_t /*k*/, std::size_t d) const {
    return d >= marked_.size() || marked_[d];
  }

  void mark(std::size_t d) {
    if (d < marked_.size())
      marked_[d] = 1;
  }

  // Postprocess: clear all marks for reuse.
  void clear() {
    sdsl::util::set_to_value(marked_, 0);
  }

 private:
  TBV marked_;
};

//~~~~~~~


// No-op preprocess hook for SADA (no coordinate transformation needed).
struct PreprocessNoOp {
  void operator()(std::size_t&, std::size_t&) const {}
};

//~~~~~~~


// Translates SA-space [sp, ep) into run-space via run_heads rank, and stores
// the original SA-space range on the ILCP state so the report fan-out can
// recover it to determine per-run SA boundaries.
//
// TIlcpState must expose:
//   setInitialRange(sp, ep) — stash the original SA-space half-open range
//   rank()                  — reference to run_heads rank structure
template <typename TIlcpState>
class PreprocessILCP {
 public:
  explicit PreprocessILCP(TIlcpState& state) : state_{state} {}

  void operator()(std::size_t& sp, std::size_t& ep) const {
    state_.setInitialRange(sp, ep);
    const auto& rank = state_.rank();
    sp = rank(sp + 1) - 1;
    ep = rank(ep);
  }

 private:
  TIlcpState& state_;
};

//~~~~~~~


// Postprocess hook: resets the MarkedReported bitvector after a query.
// Use when you want to reuse the same MarkedReported across multiple queries
// without per-query allocation.
template <typename TBV>
class PostprocessClearMarked {
 public:
  explicit PostprocessClearMarked(MarkedReported<TBV>& mr) : mr_{mr} {}

  void operator()() const {
    mr_.clear();
  }

 private:
  MarkedReported<TBV>& mr_;
};

}  // namespace dret::rmq
