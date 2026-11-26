//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Document frequencies against a direct count: for every pattern and document, the number of
// occurrences of the pattern in the document's text (overlapping ones included), over every index
// that serves frequencies (doc_list_index.h, FrequencyCore), on versions of one random reference with
// a few short documents and one-character ones, so ranges hold repeated documents, documents once, and
// patterns that cross no delimiter but end at one.
//

#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <map>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sr-index/r_index.h"

#include "dret/backend/da.h"
#include "dret/backend/rlcsa.h"
#include "dret/backend/sa_phi.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_brute_backend.h"
#include "dret/doc_list/doc_list_slp.h"

#include "base_test.h"

namespace {

using S = std::reference_wrapper<dret::GenericStorage>;
using A = dret::Alphabet<>;
using RC = sri::RIndexCount<S, A>;

template <typename TIndex>
struct Kind {
  using Index = TIndex;
  static Index Make(S t_storage) { return Index(t_storage); }
};

struct BruteRI : Kind<dret::DocListIdxBrute<S, A, sri::RIndex<S>>> {};
struct BruteSRI : Kind<dret::DocListIdxBrute<S, A>> {};
struct BruteDA : Kind<dret::DocListIdxBruteBackend<S>> {};
struct BruteSAPhi : Kind<dret::DocListIdxBruteBackend<S, A, RC, dret::backend::SAPhi<S, A::int_width>>> {};
struct BruteRLCSA : Kind<dret::DocListIdxBruteBackend<S, A, RC, dret::backend::RLCSA<S, A::int_width>>> {};
struct BruteGCDA : Kind<dret::DocListIdxSLP<S>> {};
using Kinds = testing::Types<BruteRI, BruteSRI, BruteDA, BruteSAPhi, BruteRLCSA, BruteGCDA>;

using Freqs = std::map<std::size_t, std::size_t>;

template <typename TKind>
class DocFreqTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    std::mt19937 rng(20261126);
    std::uniform_int_distribution<int> base(0, 3), pos(0, 299), coin(0, 9);
    std::string ref;
    for (int i = 0; i < 300; ++i) ref += "ACGT"[base(rng)];
    for (int d = 0; d < 24; ++d) {
      std::string doc;
      if (d % 8 == 7) {
        doc = std::string(1, "ACGT"[base(rng)]);  // one character
      } else {
        doc = ref;
        for (int k = 0; k < 6; ++k) doc[pos(rng)] = "ACGT"[base(rng)];
        if (coin(rng) < 3) doc = doc.substr(0, 50 + pos(rng) / 3);
      }
      docs_.push_back(doc);
      data_ += doc;
      data_ += '\1';
    }
    Init(data_);

    std::set<std::string> patterns{"A", "C", "G", "T"};
    std::uniform_int_distribution<std::size_t> start(0, data_.size() - 1), len(1, 10);
    while (patterns.size() < 400) {
      auto p = data_.substr(start(rng), len(rng));
      if (p.find('\1') == std::string::npos) patterns.insert(p);
    }
    patterns_.assign(patterns.begin(), patterns.end());
  }

  // The occurrences of t_pattern in each document, overlapping ones included.
  Freqs Expected(const std::string& t_pattern) const {
    Freqs f;
    for (std::size_t d = 0; d < docs_.size(); ++d) {
      std::size_t n = 0;
      for (auto at = docs_[d].find(t_pattern); at != std::string::npos; at = docs_[d].find(t_pattern, at + 1)) ++n;
      if (n) f[d] = n;
    }
    return f;
  }

  std::string data_;
  std::vector<std::string> docs_;
  std::vector<std::string> patterns_;
};

TYPED_TEST_SUITE(DocFreqTest, Kinds);

TYPED_TEST(DocFreqTest, MatchesADirectCount) {
  dret::GenericStorage storage;
  auto index = TypeParam::Make(std::ref(storage));
  construct(index, this->config_);

  for (const auto& p : this->patterns_) {
    Freqs got;
    std::size_t reports = 0;
    index.Frequencies(p, [&](std::size_t d, std::size_t f) {
      got[d] += f;
      ++reports;
    });
    EXPECT_EQ(reports, got.size()) << "a document reported twice for '" << p << "'";
    EXPECT_EQ(got, this->Expected(p)) << "pattern '" << p << "'";
  }
}

}  // namespace
