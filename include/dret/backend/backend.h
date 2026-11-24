//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// Document-array backends: what reads the document of a suffix-array position. Every index that
// expands positions goes through one -- the brute-force indexes over the whole range, GCDA and PDL
// at the ends a sampled node does not cover, the RMQ cores at the runs they visit:
//
//   PlainDA          the plain document array                        (backend/da.h)
//   SAPhi            SA[i] by Phi walks of the r-index or sr-index,
//                    then the document ends                          (backend/sa_phi.h)
//   RLCSA            the RLCSA's locate                              (backend/rlcsa.h)
//   GCDAGrammar      the sampled GCDA grammar                        (backend/grammar.h)
//   BareGrammar      a GCDA grammar without its sampled tree         (backend/grammar.h)
//   DiffGCDAGrammar  the sampled differential grammar                (backend/grammar.h)
//
// A backend is an IndexBaseWithExternalStorage (it loads its structures from the shared storage),
// has a construct(backend, config) that builds them, and answers:
//
//   operator()(i)            the document of position i
//   operator()(b, e, r)      r(d) for the document of every position of [b, e), in no given order
//   ExpandUntil(b, e, f)     f(d) for [b, e) in order, until f returns false
//
// A backend may also find the range of a pattern itself (Count, as the RLCSA does, whose suffix array
// orders suffixes differently around the document delimiters, so its positions are its own): an index
// over it then uses it in place of the r-index (CountsPattern).
//

#pragma once

#include <concepts>
#include <cstddef>
#include <utility>

namespace dret::backend {

namespace internal {
// A callable that accepts a document and returns whether to go on, for the ExpandUntil check.
struct AnyUntil {
  bool operator()(std::size_t) const { return true; }
};
// A report callable, for the range-expansion check.
struct AnyReport {
  void operator()(std::size_t) const {}
};
}  // namespace internal

template <typename T>
concept DocArrayBackend = requires(const T& t_backend, std::size_t t_i, internal::AnyReport& t_report) {
  { t_backend(t_i) } -> std::convertible_to<std::size_t>;
  t_backend(t_i, t_i, t_report);
  t_backend.ExpandUntil(t_i, t_i, internal::AnyUntil{});
};

// The backend finds the suffix-array range of a pattern itself.
template <typename T, typename TPattern>
concept CountsPattern = requires(const T& t_backend, const TPattern& t_pattern) {
  { t_backend.Count(t_pattern) } -> std::convertible_to<std::pair<std::size_t, std::size_t>>;
};

}  // namespace dret::backend
