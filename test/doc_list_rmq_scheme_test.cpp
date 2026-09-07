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
// The CILCP marker-stop regression that used to live here has moved to
// test/cilcp_fanout_test.cpp, where it runs against the real cores instead of
// against a hand-built run array. It belongs there: the defect was never in
// this scheme but in the range the caller recurses over, and the fix (restrict
// the marker stop to the runs contained in the query range, expand the two
// boundary runs separately) lives in IlcpLikeLeanCore::findDocs. A scheme-level
// test could only pin behaviour no core asks for.
//
// The shape that breaks the stop, for reference: run values [5, 5, 7, 3, 9]
// with leftmost docs [X, Y, Z, Y, W]. Y leads both run 1 and run 3, and run 3
// holds the smaller value, so it wins the RMQ over the whole range, reports Y,
// and stops the recursion before run 2 ever surrenders Z.
