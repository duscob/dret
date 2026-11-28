//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Randomised differential test of the document-frequency indexes against a direct count, on random
// collections of 3-7 documents of 1-10 symbols over {A, B, C}: the RMQ partitions then have runs that
// straddle a query range endpoint and documents a single position wide, which is where the leftmost and
// rightmost passes of the RMQ frequency cores (doc_freq/rmq_core.h) and the reversal they rest on can
// go wrong. Deterministic: DRET_DIFF_SEED and DRET_DIFF_COLLECTIONS override the seed and the count.
//
//   DRET_DIFF_COLLECTIONS=700 ./doc_freq_random_diff_test
//

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <functional>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "dret/backend/diff_sa.h"
#include "dret/backend/wavelet_da.h"
#include "dret/doc_freq/oracle.h"
#include "dret/doc_freq/rmq_core.h"
#include "dret/doc_list/doc_list_brute_backend.h"
#include "dret/doc_list/doc_list_gcda.h"

#include "base_test.h"

namespace {

using S = std::reference_wrapper<dret::GenericStorage>;
using A = dret::Alphabet<>;
using RC = sri::RIndexCount<S, A>;
constexpr auto kW = A::int_width;

using Freqs = std::map<std::size_t, std::size_t>;

Freqs DirectCount(const std::vector<std::string>& t_docs, const std::string& t_pattern) {
  Freqs f;
  for (std::size_t d = 0; d < t_docs.size(); ++d) {
    std::size_t n = 0;
    for (auto at = t_docs[d].find(t_pattern); at != std::string::npos; at = t_docs[d].find(t_pattern, at + 1)) ++n;
    if (n) f[d] = n;
  }
  return f;
}

std::vector<std::string> AllPatterns(const std::vector<std::string>& t_docs) {
  std::set<std::string> pats;
  for (const auto& doc : t_docs)
    for (std::size_t i = 0; i < doc.size(); ++i)
      for (std::size_t m = 1; m <= 4 && i + m <= doc.size(); ++m) pats.insert(doc.substr(i, m));
  return {pats.begin(), pats.end()};
}

std::vector<std::string> RandomCollection(std::mt19937& t_gen) {
  std::uniform_int_distribution<std::size_t> n_docs(3, 7), doc_len(1, 10);
  std::uniform_int_distribution<int> letter(0, 2);
  std::vector<std::string> docs(n_docs(t_gen));
  for (auto& doc : docs) {
    doc.resize(doc_len(t_gen));
    for (auto& c : doc) c = static_cast<char>('A' + letter(t_gen));
  }
  return docs;
}

std::size_t EnvOr(const char* t_name, std::size_t t_fallback) {
  if (const char* v = std::getenv(t_name)) {
    const auto n = std::strtoull(v, nullptr, 10);
    if (n > 0) return static_cast<std::size_t>(n);
  }
  return t_fallback;
}

template <typename TIndex>
struct Kind {
  using Index = TIndex;
  static Index Make(S t_storage) { return Index(t_storage); }
};

using WT = dret::backend::WaveletDA<S, kW>;
using DSA = dret::backend::DiffSA<S, kW>;
using Rank = dret::freq::WtRankOracle<S, kW>;
using Isa = dret::freq::DocIsaOracle<S, kW>;
using Cores = dret::freq::Cores<S, kW>;

template <template <typename> class TCoreOf>
struct RmqWT : Kind<dret::freq::DocFreqIdxRMQ<S, A, RC, WT, TCoreOf, Rank>> {};
// Sample spacing 4, so the oracle's lookups land between samples.
template <template <typename> class TCoreOf>
struct RmqIsa : Kind<dret::freq::DocFreqIdxRMQ<S, A, RC, DSA, TCoreOf, Isa>> {
  using Index = dret::freq::DocFreqIdxRMQ<S, A, RC, DSA, TCoreOf, Isa>;
  static Index Make(S t_storage) { return Index(t_storage, 4u); }
};
struct GCDA : Kind<dret::gcda::DocFreqIdxGCDA<S>> {
  static Index Make(S t_storage) { return Index(t_storage, 4, 2); }
};
struct BruteDA : Kind<dret::DocListIdxBruteBackend<S>> {};

using Kinds = testing::Types<RmqWT<Cores::Sada>,
                             RmqWT<Cores::Ilcp>,
                             RmqWT<Cores::Cilcp>,
                             RmqIsa<Cores::Sada>,
                             RmqIsa<Cores::Ilcp>,
                             RmqIsa<Cores::Cilcp>,
                             GCDA,
                             BruteDA>;

template <typename TKind>
class DocFreqRandomDiffTest : public BaseConfigTests<8> {
 protected:
  // One directory per collection: Init reuses the test's, so wipe it first.
  void InitFresh(const std::string& t_data) {
    if (!tmp_dir_.empty()) {
      sdsl::util::delete_all_files(config_.file_map);
      std::error_code ec;
      std::filesystem::remove_all(tmp_dir_, ec);
    }
    Init(t_data);
  }
};

TYPED_TEST_SUITE(DocFreqRandomDiffTest, Kinds);

TYPED_TEST(DocFreqRandomDiffTest, CountsEveryOccurrence) {
  const auto seed = EnvOr("DRET_DIFF_SEED", 20261125u);
  const auto n_collections = EnvOr("DRET_DIFF_COLLECTIONS", 200u);
  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  std::size_t n_queries = 0, n_wrong = 0;
  for (std::size_t c = 0; c < n_collections; ++c) {
    const auto docs = RandomCollection(gen);
    std::string data;
    for (const auto& d : docs) data += d + '\1';
    this->InitFresh(data);

    dret::GenericStorage storage;  // a fresh one: a shared storage would hand out the last collection's
    auto index = TypeParam::Make(std::ref(storage));
    construct(index, this->config_);

    for (const auto& p : AllPatterns(docs)) {
      Freqs got;
      index.Frequencies(p, [&got](std::size_t d, std::size_t f) { got[d] += f; });
      ++n_queries;
      const auto want = DirectCount(docs, p);
      if (got != want && ++n_wrong <= 5) {
        ADD_FAILURE() << "collection " << c << " (seed " << seed << "), pattern \"" << p << "\"\n"
                      << "  docs: " << testing::PrintToString(docs) << "\n"
                      << "  got:  " << testing::PrintToString(got) << "\n"
                      << "  want: " << testing::PrintToString(want);
      }
    }
  }
  EXPECT_EQ(n_wrong, 0u) << n_wrong << " wrong answers out of " << n_queries << " queries";
  std::cout << "[          ] " << n_queries << " queries over " << n_collections << " collections (seed " << seed
            << ")" << std::endl;
}

}  // namespace
