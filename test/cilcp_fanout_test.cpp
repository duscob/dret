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

// Witness A -- breaks CILCP.  Pattern "BC" has SA range [15, 19) with
// DA = [1, 0, 1, 1] and ILCP = [1, 1, 2, 3].
const std::vector<std::string> kDocsA = {"BCAB", "BCBCBBCAA", "A", "BA", "B"};

// Witness B -- breaks CILCP-S.  Pattern "BB" has SA range [11, 14) with
// DA = [1, 0, 1] and ILCP = [1, 1, 2].
const std::vector<std::string> kDocsB = {"BB", "BCCCBBCBB", "BAB", "A"};

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
                             dret::rmq::SadaCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::SadaSCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::IlcpSCore<ExternalGenericStorage>>,
    dret::rmq::DocListIdxRMQ<ExternalGenericStorage, dret::Alphabet<>,
                             sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                             dret::rmq::CilcpSCore<ExternalGenericStorage>>>;

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
