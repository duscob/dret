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
//
// kStopOnReported gates the marker-based recursion-stop. With true (default),
// recursion halts whenever a subrange's RMQ-min has an already-reported doc —
// a real ~10x speedup that is empirically sound for SADA (prev_doc-based RMQ)
// and ILCP (ilcp-constant runs). With false, recursion always explores both
// halves; needed for CILCP, whose Rule-2 RLE merges by doc-equality with min
// tracking, breaking the marker-based stop invariant and silently dropping
// docs in subranges whose RMQ-min coincides with an earlier-marked doc.
//
// kAlwaysReport gates whether report() is invoked even when the RMQ-min's doc
// is already reported. For SADA one RMQ position carries exactly one document,
// so gating on the marker loses nothing. For the ILCP family a position stands
// for a whole RUN, and report() is also what fans the run out; a run whose head
// doc happens to be a duplicate can still hold the FIRST occurrence of other
// documents, and gating drops them. That is sound only when the run's own
// values pin down its contents:
//   - ILCP / ILCP-S: runs are ilcp-constant, so a run reachable past the stop
//     test consists purely of first occurrences and its head doc cannot already
//     be reported — gating is a no-op.
//   - CILCP / CILCP-S: runs merge by document, so a run's stored value is a
//     MIN and may even be attained outside [sp, ep). Such a run can report a
//     duplicate occurrence early and thereby suppress the fan-out of a later
//     run that holds a genuinely new document. These cores must always report.
// report() is responsible for de-duplicating its own head document.
template <bool kStopOnReported = true,
          bool kAlwaysReport = false,
          typename TRMQ, typename TGetDoc, typename TIsReported, typename TReport>
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
  const bool already = is_reported(k, d);
  if (!already || kAlwaysReport) {
    report(k, d);
  }
  if constexpr (kStopOnReported) {
    if (already) return;
  }
  ListDocsRMQScheme<kStopOnReported, kAlwaysReport>(bp, k, rmq, get_doc, is_reported, report);
  ListDocsRMQScheme<kStopOnReported, kAlwaysReport>(k + 1, ep, rmq, get_doc, is_reported, report);
}

// Canonical Sadakane-style LEFTMOST RMQ traversal: recursion stops by an
// explicit value-based predicate, not by marker state. This is the algorithm
// used by the SADA, ILCP and CILCP★ papers (Sadakane 2007; Gagie, Navarro,
// Puglisi 2014; Cobas, Mäkinen, Rossi SPIRE 2020).
//
// stop_pred(k, bp) -> bool returns true when the RMQ-min over the current
// subrange [bp, ep) cannot contribute any new leftmost-doc occurrence:
//   - SADA-S:   prev_doc[k] >= bp_subrange (no doc in [bp..k] has its
//               previous occurrence before bp).
//   - ILCP-S / CILCP-S: run_values[k] >= m (no position in the subrange
//               has ILCP < m, so Lemma 1 / Lemma 2 give no leftmost-doc
//               here).
//
// is_reported still gates the *emission* of doc d (same MarkedReported as
// the original scheme) — multiple runs can share a leftmost-doc, and we
// don't want to double-emit. It no longer gates the recursion.
// kAlwaysReport carries the same meaning as in ListDocsRMQScheme above: the
// ILCP-family cores whose runs merge by document (CILCP-S) must fan a visited
// run out even when its head document is already reported.
template <bool kAlwaysReport = false,
          typename TRMQ, typename TGetDoc, typename TStopPred,
          typename TIsReported, typename TReport>
void ListDocsRMQSchemeDepth(std::size_t bp,
                            std::size_t ep,
                            const TRMQ& rmq,
                            const TGetDoc& get_doc,
                            const TStopPred& stop_pred,
                            const TIsReported& is_reported,
                            TReport& report) {
  if (bp >= ep)
    return;
  const auto k = rmq(bp, ep - 1);
  if (stop_pred(k, bp)) return;
  const auto d = get_doc(k);
  if (!is_reported(k, d) || kAlwaysReport) {
    report(k, d);
  }
  ListDocsRMQSchemeDepth<kAlwaysReport>(bp, k, rmq, get_doc, stop_pred, is_reported, report);
  ListDocsRMQSchemeDepth<kAlwaysReport>(k + 1, ep, rmq, get_doc, stop_pred, is_reported, report);
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
