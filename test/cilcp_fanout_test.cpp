//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Differential test for the ILCP-family run fan-out. The collection below is a
// minimal witness: pattern "BC" has SA range [15, 19) with DA = [1, 0, 1, 1] and
// ILCP = [1, 1, 2, 3]. The CILCP run [17..19] is a same-doc run for doc 1 whose
// stored minimum (0) is attained at position 19, OUTSIDE the query range, so it
// wins the RMQ and reports doc 1 from a duplicate occurrence. The run [15..16]
// (DA = [1, 0]) then has an already-reported head doc, so its fan-out never
// runs and doc 0 is never listed.

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_rmq.h"

#include "base_test.h"

using ExternalGenericStorage = sri::GenericStorage;

namespace {

class Collector : public std::vector<std::size_t> {
 public:
  void operator()(std::size_t d) { emplace_back(d); }
  void operator()() {
    std::sort(begin(), end());
    erase(std::unique(begin(), end()), end());
  }
};

std::vector<std::size_t> BruteForce(const std::vector<std::string>& docs, const std::string& p) {
  std::set<std::size_t> s;
  for (std::size_t d = 0; d < docs.size(); ++d)
    if (docs[d].find(p) != std::string::npos) s.insert(d);
  return {s.begin(), s.end()};
}

std::vector<std::string> AllPatterns(const std::vector<std::string>& docs) {
  std::set<std::string> pats;
  for (const auto& doc : docs)
    for (std::size_t i = 0; i < doc.size(); ++i)
      for (std::size_t m = 1; m <= 4 && i + m <= doc.size(); ++m)
        pats.insert(doc.substr(i, m));
  return {pats.begin(), pats.end()};
}

}  // namespace

// Witness A -- breaks CILCP-L (the marker-based core).  Pattern "BC" has SA
// range [15, 19) with DA = [1, 0, 1, 1] and ILCP = [1, 1, 2, 3].
const std::vector<std::string> kDocsA = {"BCAB", "BCBCBBCAA", "A", "BA", "B"};

// Witness B -- breaks BOTH merged cores.  Pattern "BB" has SA range [11, 14)
// with DA = [1, 0, 1] and ILCP = [1, 1, 2]. It was added for CILCP, but it also
// fails a CILCP-L that recurses over the whole range, so it is the one witness
// here that still guards both.
const std::vector<std::string> kDocsB = {"BB", "BCCCBBCBB", "BAB", "A"};

// Witness C -- breaks CILCP-L, and replaces witness A in that role.
//
// Verified by deliberately breaking each core and checking which witnesses
// notice; do that again before trusting any of them, and mind the type index
// (DocListIdxBrute is 0 here, so CILCP-L is 3 and CILCP is 6).
//
// Witness A no longer bites. It was found against the CILCP partition as it
// stood on 2026-08-29; f3720ff then rebuilt that partition on the CMR20 rule,
// and under the new runs witness A passes even against a CILCP-L that recurses
// over the whole range with the fan-out gate on -- the exact defect it was
// added to catch. It is kept as a correctness case, but it guards nothing.
//
// This one was found by cilcp_random_diff_test (seed 20260907, collection 48)
// and is verified to fail against that same broken build, so it guards the
// boundary-run handling in IlcpLikeLeanCore::findDocs. Pattern "BAB" returns
// {0} where the answer is {0, 2}.
const std::vector<std::string> kDocsC = {"BBBABABAA", "A", "CBCBABAAC"};

std::string Concat(const std::vector<std::string>& docs) {
  std::string s;
  for (const auto& d : docs) { s += d; s += '\1'; }
  return s;
}

template <typename T>
class CilcpFanoutTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {}
  sri::GenericStorage storage_;
};

using CilcpFanoutTypes = ::testing::Types<
    dret::DocListIdxBrute<ExternalGenericStorage>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaLCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpLCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpLCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpCore<ExternalGenericStorage>>>;

TYPED_TEST_SUITE(CilcpFanoutTypedTests, CilcpFanoutTypes);

TYPED_TEST(CilcpFanoutTypedTests, lists_every_document_witness_a) {
  this->Init(Concat(kDocsA));
  TypeParam index(std::ref(this->storage_));
  construct(index, this->config_);

  for (const auto& p : AllPatterns(kDocsA)) {
    Collector got;
    index.Search(p, std::ref(got));
    got();
    EXPECT_THAT(got, testing::ElementsAreArray(BruteForce(kDocsA, p))) << "pattern \"" << p << '"';
  }
}

TYPED_TEST(CilcpFanoutTypedTests, lists_every_document_witness_c) {
  this->Init(Concat(kDocsC));
  TypeParam index(std::ref(this->storage_));
  construct(index, this->config_);

  for (const auto& p : AllPatterns(kDocsC)) {
    Collector got;
    index.Search(p, std::ref(got));
    got();
    EXPECT_THAT(got, testing::ElementsAreArray(BruteForce(kDocsC, p))) << "pattern \"" << p << '"';
  }
}

TYPED_TEST(CilcpFanoutTypedTests, lists_every_document_witness_b) {
  this->Init(Concat(kDocsB));
  TypeParam index(std::ref(this->storage_));
  construct(index, this->config_);

  for (const auto& p : AllPatterns(kDocsB)) {
    Collector got;
    index.Search(p, std::ref(got));
    got();
    EXPECT_THAT(got, testing::ElementsAreArray(BruteForce(kDocsB, p))) << "pattern \"" << p << '"';
  }
}

// The paper's claim is that CILCP-L and CILCP partition IDENTICALLY -- they are
// the same CMR20 runs, and differ only in whether the run values are stored.
// Enforce it here: the run_heads structure must match byte for byte, while
// CILCP must be strictly larger overall by its run_values.
class CilcpPartitionTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override { this->Init(std::string("BCAB\1BCBCBBCAA\1A\1BA\1B\1")); }
  sri::GenericStorage storage_;
};

TEST_F(CilcpPartitionTest, cilcp_and_cilcp_s_share_one_partition) {
  using Cilcp  = dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                     sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                     dret::rmq::CilcpLCore<ExternalGenericStorage>>;
  using CilcpS = dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                     sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                     dret::rmq::CilcpCore<ExternalGenericStorage>>;
  auto bytes = [](const auto& idx, const std::string& key) -> std::size_t {
    for (const auto& [k, v] : idx.core().GetSizeReport())
      if (k == key) return v;
    return 0;
  };
  Cilcp  a(std::ref(this->storage_)); construct(a, this->config_);
  CilcpS b(std::ref(this->storage_)); construct(b, this->config_);

  EXPECT_GT(bytes(a, "run_heads"), 0u);
  EXPECT_EQ(bytes(a, "run_heads"), bytes(b, "run_heads")) << "partitions differ";
  EXPECT_EQ(bytes(a, "rmq"), bytes(b, "rmq")) << "run counts differ";
  EXPECT_EQ(bytes(a, "run_values"), 0u) << "CILCP-L must not store run values";
  EXPECT_GT(bytes(b, "run_values"), 0u) << "CILCP must store run values";
}
