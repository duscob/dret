//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Stateless LEFTMOST-RMQ traversal and supporting hook classes for the
// document-listing RMQ indexes (SADA / ILCP / CILCP). Kept free of grammar /
// SLP dependencies so it can be included cheaply in any listing core.
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
// For SADA: is_reported and report read and fill the query's sink (dret/sink.h):
//   sink.Contains(d), and sink.Position(d).
// For ILCP: report also fans out the rest of the run; see IlcpLikeCore.
//
// The recursion halts whenever a subrange's RMQ-min has an already-reported
// doc — a real ~10x speedup. Its soundness is a statement about the RUNS the
// caller hands us, not about this function: the stop is licensed exactly when
// an already-reported head implies the run's stored value is at least the
// pattern length, so that the whole subrange is free of first occurrences.
// That holds for SADA (one position per RMQ entry) and for ILCP (runs are
// ilcp-constant, so a run's value cannot be borrowed from outside the query
// range). It holds for a document-aware merge only at runs CONTAINED in the
// range; a run straddling an endpoint can hold a minimum attained outside it,
// win the RMQ, report a repeat occurrence and stop the recursion with pending
// documents still in the subrange. CILCP therefore does not weaken this
// stop — it restricts the range it recurses over to the contained runs and
// expands the two boundary runs itself. See IlcpLikeCore::Query.
//
// The marker test IS the stop, so a run that fails it is neither reported nor
// descended into. For SADA that loses nothing: one RMQ position carries one
// document. For the ILCP family a position stands for a whole RUN and report()
// is also what fans the run out, so skipping it would drop any first occurrence
// deeper in the run -- but a run cannot both fail the test and hold one. An
// already-reported head forces the run's stored value to at least the pattern
// length, and every position of such a run is then a repeat. That argument needs
// the value to be a minimum over positions of the query range, which is why
// CILCP hands us only the contained runs.
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
  if (is_reported(k, d)) return;
  report(k, d);
  ListDocsRMQScheme(bp, k, rmq, get_doc, is_reported, report);
  ListDocsRMQScheme(k + 1, ep, rmq, get_doc, is_reported, report);
}

// ListDocsRMQScheme with the head lookup, the marker test and the report fused
// into visit(k): it reads the run at k, reports it unless its head document is
// already reported, and returns whether it did. The ILCP family visits a run in
// one in-order backend expansion (VisitRun in doc_list_rmq.h) instead of a head
// lookup followed by a fan-out that locates its start again. The stop is the
// marker stop above, with the same soundness conditions.
template <typename TRMQ, typename TVisit>
void ListDocsRMQSchemeVisit(std::size_t bp, std::size_t ep, const TRMQ& rmq, const TVisit& visit) {
  if (bp >= ep)
    return;
  const auto k = rmq(bp, ep - 1);
  if (!visit(k)) return;
  ListDocsRMQSchemeVisit(bp, k, rmq, visit);
  ListDocsRMQSchemeVisit(k + 1, ep, rmq, visit);
}

// LEFTMOST RMQ traversal whose recursion stops by an explicit value-based
// predicate, not by marker state: Muthukrishnan's rule (2002), with the array the
// RMQ is built over stored. dret's -V variants (SADA-V, ILCP-V, CILCP-V) use it.
// The published SADA, ILCP and CILCP (Sadakane 2007; Gagie et al. 2017; Cobas,
// Makinen and Rossi, SPIRE 2020) stop by marking instead and store no values:
// ListDocsRMQScheme above.
//
// stop_pred(k) -> bool returns true when the RMQ-min over the current subrange
// cannot contribute any new leftmost-doc occurrence. Every predicate is a
// question about k alone, against the ORIGINAL query range -- none of them
// needs the current subrange start, which is why none is passed:
//   - SADA-V:   prev_doc[k] >= sp (no position in the subrange has its previous
//               occurrence before sp, so none is a first occurrence).
//   - ILCP-V / CILCP-V: run_values[k] >= m (no position in the subrange has
//               ILCP < m, so Lemma 1 / Lemma 2 give no leftmost-doc here).
// Testing SADA-V against the subrange start instead would still be sound but
// prunes less, and would let the traversal report a position that is first in
// its subrange without being first in the range -- that is, a duplicate.
//
// Unlike the scheme above, this traversal does not consult reporting state at
// all. Only stop_pred decides where it goes, and every run it reaches is
// reported. That is not an optimization, it is what makes the scheme correct
// for ANY partition into runs: pruning is a statement about stored values, and
// a marker can be set by a run that reported a REPEAT occurrence (a run whose
// stored minimum was attained outside [sp, ep) can do exactly that), so reading
// it here would prune on evidence the values never supported.
//
// De-duplicating the emission is therefore the caller's business, and every
// core's report() does it. Measured over the differential corpus: an
// is_reported gate here would have fired 0 times in 31059 reported runs for
// ILCP-V, and 166 times in 31229 for CILCP-V, of which 5 carried a document that
// nothing else would have listed. SADA-V would have fired 15185 times in 55480,
// all of them redundant emissions its report() now absorbs.
template <typename TRMQ, typename TGetDoc, typename TStopPred, typename TReport>
void ListDocsRMQSchemeDepth(std::size_t bp,
                            std::size_t ep,
                            const TRMQ& rmq,
                            const TGetDoc& get_doc,
                            const TStopPred& stop_pred,
                            TReport& report) {
  if (bp >= ep)
    return;
  const auto k = rmq(bp, ep - 1);
  if (stop_pred(k)) return;
  report(k, get_doc(k));
  ListDocsRMQSchemeDepth(bp, k, rmq, get_doc, stop_pred, report);
  ListDocsRMQSchemeDepth(k + 1, ep, rmq, get_doc, stop_pred, report);
}

// ListDocsRMQSchemeDepth with the lookup and the report fused into visit(k), as
// in ListDocsRMQSchemeVisit; the value predicate alone still decides the descent.
template <typename TRMQ, typename TStopPred, typename TVisit>
void ListDocsRMQSchemeDepthVisit(std::size_t bp,
                                 std::size_t ep,
                                 const TRMQ& rmq,
                                 const TStopPred& stop_pred,
                                 const TVisit& visit) {
  if (bp >= ep)
    return;
  const auto k = rmq(bp, ep - 1);
  if (stop_pred(k)) return;
  visit(k);
  ListDocsRMQSchemeDepthVisit(bp, k, rmq, stop_pred, visit);
  ListDocsRMQSchemeDepthVisit(k + 1, ep, rmq, stop_pred, visit);
}

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

}  // namespace dret::rmq
