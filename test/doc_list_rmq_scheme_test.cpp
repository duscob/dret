//
// Unit tests for dret::rmq::ListDocsRMQScheme + MarkedReported.
//
// These exercise the stateless RMQ-doc-listing scheme directly with synthetic
// inputs — no index construction, no Config, no Factory. Split out of
// doc_list_test.cpp so the heavy integration suite stays focused on the full
// doc-list family.
//

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_succinct_sct.hpp>

#include "dret/rmq/doc_list_rmq_scheme.h"

//~~~~~~~
// Basic scheme coverage on a synthetic DA + SADA-style prev_doc RMQ.


class ListDocsRMQSchemeTest : public ::testing::Test {
 protected:
  // DA = {0, 1, 2, 0, 1}  (3 distinct docs, 5 SA positions)
  // prev_doc[i] = last_occ[da[i]] before position i (0 when first occurrence).
  // prev_doc = {0, 0, 0, 0, 1}
  void SetUp() override {
    const std::vector<std::size_t> da_vec = {0, 1, 2, 0, 1};
    n_doc_ = 3;

    sdsl::int_vector<> prev_doc(da_vec.size(), 0, sdsl::bits::hi(da_vec.size()) + 1);
    std::vector<std::size_t> last_occ(n_doc_, 0);
    for (std::size_t i = 0; i < da_vec.size(); ++i) {
      prev_doc[i] = last_occ[da_vec[i]];
      last_occ[da_vec[i]] = i;
    }

    da_.resize(da_vec.size());
    for (std::size_t i = 0; i < da_vec.size(); ++i)
      da_[i] = da_vec[i];

    rmq_ = sdsl::rmq_succinct_sct<true>(&prev_doc);
  }

  std::vector<std::size_t> query(std::size_t sp, std::size_t ep_closed) {
    dret::rmq::MarkedReported mr(n_doc_);
    std::vector<std::size_t> reported;
    auto get_doc = [this](std::size_t k) {
      return static_cast<std::size_t>(da_[k]);
    };
    auto report = [&mr, &reported](std::size_t /*k*/, std::size_t d) {
      mr.mark(d);
      reported.push_back(d);
    };
    dret::rmq::ListDocsRMQScheme(sp, ep_closed + 1, rmq_, get_doc, mr, report);
    std::sort(reported.begin(), reported.end());
    return reported;
  }

  sdsl::int_vector<> da_;
  sdsl::rmq_succinct_sct<true> rmq_;
  std::size_t n_doc_ = 0;
};

TEST_F(ListDocsRMQSchemeTest, full_range_reports_all_docs) {
  EXPECT_EQ(query(0, 4), (std::vector<std::size_t>{0, 1, 2}));
}

TEST_F(ListDocsRMQSchemeTest, suffix_range_reports_two_docs) {
  EXPECT_EQ(query(3, 4), (std::vector<std::size_t>{0, 1}));
}

TEST_F(ListDocsRMQSchemeTest, single_position_reports_one_doc) {
  EXPECT_EQ(query(4, 4), (std::vector<std::size_t>{1}));
}

TEST_F(ListDocsRMQSchemeTest, middle_range_deduplicates) {
  // DA[1]=1, DA[2]=2, DA[3]=0
  EXPECT_EQ(query(1, 3), (std::vector<std::size_t>{0, 1, 2}));
}

TEST_F(ListDocsRMQSchemeTest, empty_range_reports_nothing) {
  EXPECT_EQ(query(2, 1), (std::vector<std::size_t>{}));
}

//~~~~~~~
// Regression test for the CILCP-specific correctness bug discovered on the
// `page` collection: the marker-based recursion-stop (kStopOnReported=true)
// is unsound for CILCP because Rule-2 RLE merges runs by doc-equality, so a
// subrange's RMQ-min can land on an already-marked doc while pending docs
// remain in the same subrange. The kStopOnReported=false dispatch (which
// IlcpLikeCore uses for CILCP) must report every doc in the range.
//
// Setup (hand-trace from docs/cilcp_query_perf.md):
//   docs = [X, X, Y, Z, Y, W]  =  [0, 0, 1, 2, 1, 3]
//   ilcp = [10, 5, 5, 7, 3, 9]
// CILCP RLE produces 5 runs with leftmost-docs {X, Y, Z, Y, W} and
//   run_values = [5, 5, 7, 3, 9].
// Querying the full run-space [0, 5) must report {X, Y, Z, W}.
// With kStopOnReported=true the recursion stops prematurely and drops Z;
// with kStopOnReported=false (the CILCP dispatch) all four are reported.

class ListDocsRMQSchemeCilcpRegressionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // run_values mirror CILCP's per-run min(ilcp); leftmost-doc per run is
    // what get_doc(run_index) must return. Y appears as the leftmost-doc of
    // both R1 and R3, which is the structural pattern that breaks the
    // marker-based stop on the standard scheme.
    run_values_ = {5, 5, 7, 3, 9};
    run_leftmost_doc_ = {0, 1, 2, 1, 3};  // X, Y, Z, Y, W
    n_doc_ = 4;
    rmq_ = sdsl::rmq_succinct_sct<true>(&run_values_);
  }

  template <bool kStopOnReported>
  std::vector<std::size_t> query(std::size_t bp, std::size_t ep) {
    dret::rmq::MarkedReported mr(n_doc_);
    std::vector<std::size_t> reported;
    auto get_doc = [this](std::size_t k) { return run_leftmost_doc_[k]; };
    auto report = [&mr, &reported](std::size_t /*k*/, std::size_t d) {
      mr.mark(d);
      reported.push_back(d);
    };
    dret::rmq::ListDocsRMQScheme<kStopOnReported>(bp, ep, rmq_, get_doc, mr, report);
    std::sort(reported.begin(), reported.end());
    return reported;
  }

  std::vector<std::size_t> run_values_;
  std::vector<std::size_t> run_leftmost_doc_;
  sdsl::rmq_succinct_sct<true> rmq_;
  std::size_t n_doc_ = 0;
};

TEST_F(ListDocsRMQSchemeCilcpRegressionTest, cilcp_dispatch_reports_all_docs) {
  // The CILCP dispatch (kStopOnReported=false) must enumerate every distinct
  // leftmost-doc in the run-space range.
  EXPECT_EQ(query<false>(0, 5), (std::vector<std::size_t>{0, 1, 2, 3}));
}

TEST_F(ListDocsRMQSchemeCilcpRegressionTest, default_dispatch_drops_doc_on_this_input) {
  // Documents the bug: with the default marker-based stop, this exact input
  // misses doc Z (=2). The test asserts the observed wrong behaviour so that
  // any future change to the default scheme's stop condition surfaces here
  // and prompts a re-evaluation of the IlcpLikeCore dispatch (currently set
  // to kStopOnReported = (kVariant != CILCP)).
  EXPECT_EQ(query<true>(0, 5), (std::vector<std::size_t>{0, 1, 3}));
}
