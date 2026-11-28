//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// The listing sink (the contract is dret/sink.h): each document reported once, the first time it is
// seen. DocIndex::Search queries through it.
//

#pragma once

#include <cstddef>

#include "dret/doc_list/distinct_docs.h"
#include "dret/sink.h"

namespace dret {

// Document listing: reports each document once, to t_report, the first time it is seen, from a
// position or a list. Contains(d) tells whether d was reported, for the marking stop of the RMQ cores.
// t_n_doc sizes the bitmap; identifiers past it grow it (distinct_docs.h), so 0 is valid when the
// number of documents is not at hand.
template <typename TReport>
class ListSink {
 public:
  ListSink(std::size_t t_n_doc, const TReport& t_report) : seen_(t_n_doc), report_(&t_report) {}

  // Inlined into the cores' inner loops (an RMQ run's fan-out reports every position): a call per
  // position measured 6% on ILCP over the DA (revision-small).
  [[gnu::always_inline]] void Position(std::size_t t_doc) {
    if (seen_.Insert(t_doc)) {
#ifdef DRET_DOC_LIST_PROFILE
      ++reported_;
#endif
      (*report_)(t_doc);
    }
  }

  [[gnu::always_inline]] void Listed(std::size_t t_doc) { Position(t_doc); }

  [[gnu::always_inline]] void operator()(std::size_t t_doc) { Position(t_doc); }

  static constexpr bool Done() { return false; }

  bool Contains(std::size_t t_doc) const { return seen_.Contains(t_doc); }

#ifdef DRET_DOC_LIST_PROFILE
  // Documents reported, for the query profile (dret/doc_list/search_profile.h).
  std::size_t Reported() const { return reported_; }
#endif

 private:
  DistinctDocs seen_;
  const TReport* report_;  // a pointer, not a reference: see distinct_docs.h on register packing
#ifdef DRET_DOC_LIST_PROFILE
  std::size_t reported_ = 0;
#endif
};

}  // namespace dret
