//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Randomised differential test for the whole RMQ document-listing family.
// Every core is checked against brute force on random collections over a tiny
// alphabet, which is what makes the ILCP-family run structure interesting: short
// documents over {A, B, C} produce many multi-document runs and, crucially, many
// runs that straddle a query range endpoint. Those boundary runs are where a
// merged partition's stored minimum can be attained outside the range, and they
// are what the hand-picked witnesses in cilcp_fanout_test.cpp isolate.
//
// This is the harness to run before changing any recursion-stop or fan-out
// gate in doc_list_rmq.h. It is deterministic: the seed and the collection
// count are fixed unless overridden by DRET_DIFF_SEED / DRET_DIFF_COLLECTIONS,
// so a failure is always reproducible.
//
//   DRET_DIFF_COLLECTIONS=700 ./cilcp_random_diff_test
//

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>
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

std::string Concat(const std::vector<std::string>& docs) {
  std::string s;
  for (const auto& d : docs) { s += d; s += '\1'; }
  return s;
}

std::size_t EnvOr(const char* name, std::size_t fallback) {
  if (const char* v = std::getenv(name)) {
    const auto n = std::strtoull(v, nullptr, 10);
    if (n > 0) return static_cast<std::size_t>(n);
  }
  return fallback;
}

// A three-letter alphabet keeps occ/ndoc high and produces the long same-document
// runs the merge is built to collapse. Documents as short as one character matter
// too: both hand-picked witnesses need them, because a one-character document
// puts a run endpoint inside a tiny suffix-array range and that is how a run comes
// to straddle a query boundary.
std::vector<std::string> RandomCollection(std::mt19937& gen) {
  std::uniform_int_distribution<std::size_t> n_docs(3, 7);
  std::uniform_int_distribution<std::size_t> doc_len(1, 10);
  std::uniform_int_distribution<int> letter(0, 2);

  std::vector<std::string> docs(n_docs(gen));
  for (auto& doc : docs) {
    doc.resize(doc_len(gen));
    for (auto& c : doc) c = static_cast<char>('A' + letter(gen));
  }
  return docs;
}

}  // namespace

template <typename T>
class RandomDiffTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {}

  // BaseConfigTests::Init reuses one directory per test, so wipe it between
  // collections; otherwise the second collection loads the first one's caches.
  void InitFresh(const std::string& data) {
    if (!tmp_dir_.empty()) {
      sdsl::util::delete_all_files(config_.file_map);
      std::error_code ec;
      std::filesystem::remove_all(tmp_dir_, ec);
    }
    this->Init(data);
  }

  sri::GenericStorage storage_;
};

using RandomDiffTypes = ::testing::Types<
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

TYPED_TEST_SUITE(RandomDiffTypedTests, RandomDiffTypes);

TYPED_TEST(RandomDiffTypedTests, lists_every_document) {
  const auto seed = EnvOr("DRET_DIFF_SEED", 20260907u);
  const auto n_collections = EnvOr("DRET_DIFF_COLLECTIONS", 200u);

  // Seeded per type, so every core sees the same collections in the same order
  // and a failure can be compared across cores.
  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  std::size_t n_queries = 0;
  std::size_t n_wrong = 0;
  for (std::size_t c = 0; c < n_collections; ++c) {
    const auto docs = RandomCollection(gen);
    this->InitFresh(Concat(docs));

    TypeParam index(std::ref(this->storage_));
    construct(index, this->config_);

    for (const auto& p : AllPatterns(docs)) {
      Collector got;
      index.Search(p, std::ref(got));
      got();
      ++n_queries;
      const auto want = BruteForce(docs, p);
      if (got != want) {
        ++n_wrong;
        // Report the first few in full: the collection is the reproduction.
        if (n_wrong <= 5) {
          ADD_FAILURE() << "collection " << c << " (seed " << seed << "), pattern \"" << p << "\"\n"
                        << "  docs: " << testing::PrintToString(docs) << '\n'
                        << "  got:  " << testing::PrintToString(got) << '\n'
                        << "  want: " << testing::PrintToString(want);
        }
      }
    }
  }
  EXPECT_EQ(n_wrong, 0u) << n_wrong << " wrong answers out of " << n_queries << " queries over "
                         << n_collections << " collections (seed " << seed << ")";
  // Print the scale: a differential test that silently shrank its own input
  // would look exactly like a passing one.
  std::cout << "[          ] " << n_queries << " queries over " << n_collections
            << " collections (seed " << seed << ")" << std::endl;
}
