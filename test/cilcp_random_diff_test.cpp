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
#include <cstdint>
#include <type_traits>
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
#include "dret/rmq/rmq_get_doc_policies.h"
#include "dret/rmq/rmq_get_doc_sa_phi.h"

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

using S = ExternalGenericStorage;
using TCount = sri::RIndexCount<S, dret::Alphabet<>>;
constexpr auto kW = dret::Alphabet<>::int_width;

// The grammar backends answer a run head that falls inside a sampled leaf (or
// between differential samples) differently from one at its start, so they are
// built at a tiny block size: even these short collections then span many
// leaves and samples, and the mid-leaf cases are the common ones.
template <typename TCore, std::uint32_t kBlock = 4>
class WithBlock : public dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, TCount, TCore> {
  using Base = dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, TCount, TCore>;

 public:
  explicit WithBlock(const S& t_storage) : Base(t_storage, TCore(t_storage, kBlock, 2.0f)) {}
};

template <typename TCore>
using Plain = dret::rmq::DocListIdxRMQ<S, dret::Alphabet<>, TCount, TCore>;

using TBareDiff = dret::DifferentialSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>;
using GDCached = dret::rmq::GetDocSLP<S, kW>;
using GDDiffTree = dret::rmq::GetDocDSLP<S, kW>;
using GDDiff = dret::rmq::GetDocSLP_NS<S, kW, TBareDiff>;
using GDSAPhi = dret::rmq::GetDocSAPhi<S, kW>;

template <typename TGetDoc>
using IlcpL = dret::rmq::IlcpLCore<S, kW, sdsl::sd_vector<>, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
template <typename TGetDoc>
using CilcpL = dret::rmq::CilcpLCore<S, kW, sdsl::sd_vector<>, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
template <typename TGetDoc>
using Ilcp = dret::rmq::IlcpCore<S, kW, sdsl::sd_vector<>, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
template <typename TGetDoc>
using Cilcp = dret::rmq::CilcpCore<S, kW, sdsl::sd_vector<>, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;

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
                             dret::rmq::CilcpCore<ExternalGenericStorage>>,
    // The run fan-out over every other backend: the ILCP family is what expands
    // runs, so it is what a backend's range expansion can get wrong.
    WithBlock<IlcpL<GDCached>>, WithBlock<CilcpL<GDCached>>, WithBlock<Ilcp<GDCached>>, WithBlock<Cilcp<GDCached>>,
    WithBlock<IlcpL<GDDiffTree>>, WithBlock<CilcpL<GDDiffTree>>, WithBlock<Ilcp<GDDiffTree>>, WithBlock<Cilcp<GDDiffTree>>,
    WithBlock<IlcpL<GDDiff>>, WithBlock<CilcpL<GDDiff>>, WithBlock<Ilcp<GDDiff>>, WithBlock<Cilcp<GDDiff>>,
    Plain<IlcpL<GDSAPhi>>, Plain<CilcpL<GDSAPhi>>, Plain<Ilcp<GDSAPhi>>, Plain<Cilcp<GDSAPhi>>>;

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

//~~~~~~~  ExpandUntil: the in-order, early-stopping contract of every backend  ~~~~~~~
//
// The listing test above cannot see the ORDER in which a backend emits a run:
// an ILCP-family run is either all first occurrences or all repeats, so the
// documents listed do not depend on it. The order is what lets a visit stop at
// the run's head without expanding the rest (VisitRun), so it is checked here
// directly: for every range [b, e) of the document array and every stop point,
// ExpandUntil must emit exactly DA[b], DA[b + 1], ... and nothing after the stop.

// A backend plus an ILCP-L index over it, built first so the caches it loads
// from exist.
template <typename TGetDoc, typename TIndex>
struct Backend {
  using GetDoc = TGetDoc;
  using Index = TIndex;
  static TGetDoc Make(const S& t_storage) {
    if constexpr (std::is_constructible_v<TGetDoc, const S&, std::uint32_t, float>)
      return TGetDoc(t_storage, 4, 2.0f);
    else
      return TGetDoc(t_storage);
  }
};

template <typename T>
class ExpandUntilTypedTests : public RandomDiffTypedTests<T> {};

using ExpandUntilTypes = ::testing::Types<
    Backend<dret::rmq::GetDocDA<S, kW>, Plain<IlcpL<dret::rmq::GetDocDA<S, kW>>>>,
    Backend<GDCached, WithBlock<IlcpL<GDCached>>>,
    Backend<GDDiffTree, WithBlock<IlcpL<GDDiffTree>>>,
    Backend<GDDiff, WithBlock<IlcpL<GDDiff>>>,
    Backend<GDSAPhi, Plain<IlcpL<GDSAPhi>>>>;

TYPED_TEST_SUITE(ExpandUntilTypedTests, ExpandUntilTypes);

TYPED_TEST(ExpandUntilTypedTests, emits_in_order_and_stops) {
  const auto seed = EnvOr("DRET_DIFF_SEED", 20260907u);
  const auto n_collections = EnvOr("DRET_DIFF_COLLECTIONS", 200u) / 4;
  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  std::size_t n_checks = 0;
  std::size_t n_wrong = 0;
  for (std::size_t c = 0; c < n_collections; ++c) {
    const auto docs = RandomCollection(gen);
    this->InitFresh(Concat(docs));

    typename TypeParam::Index index(std::ref(this->storage_));
    construct(index, this->config_);

    dret::rmq::GetDocDA<S, kW> truth(std::ref(this->storage_));
    truth.load(this->config_);
    auto get_doc = TypeParam::Make(std::ref(this->storage_));
    get_doc.load(this->config_);

    const std::size_t n = truth.size();
    for (std::size_t b = 0; b < n; ++b) {
      for (std::size_t e = b + 1; e <= n; ++e) {
        for (std::size_t stop : {std::size_t{1}, std::size_t{2}, e - b}) {
          std::vector<std::size_t> got;
          get_doc.ExpandUntil(b, e, [&got, stop](std::size_t d) {
            got.push_back(d);
            return got.size() < stop;
          });
          std::vector<std::size_t> want;
          for (std::size_t i = b; i < e && want.size() < stop; ++i)
            want.push_back(truth(i));
          ++n_checks;
          if (got != want && ++n_wrong <= 5) {
            ADD_FAILURE() << "collection " << c << ", [" << b << ", " << e << "), stop " << stop << '\n'
                          << "  docs: " << testing::PrintToString(docs) << '\n'
                          << "  got:  " << testing::PrintToString(got) << '\n'
                          << "  want: " << testing::PrintToString(want);
          }
        }
      }
    }
  }
  EXPECT_EQ(n_wrong, 0u) << n_wrong << " wrong expansions out of " << n_checks;
  std::cout << "[          ] " << n_checks << " expansions over " << n_collections << " collections" << std::endl;
}
