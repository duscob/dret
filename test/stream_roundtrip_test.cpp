//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Every listing family through its serialization: build an index, serialize it to a stream, load a
// fresh one from that stream alone (a new storage, the cache directory deleted), and require the same
// answers, the same size report and the same size, with the stream read to its end, and the same
// bytes when the loaded index is serialized again (no path or other trace of the cache in them). The load reads
// the structures back in the order serialize() wrote them, which is what a component's field list
// (dret/fields.h) guarantees and what hand-written methods used to keep in step by hand.
//

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "sr-index/r_index.h"

#include "dret/backend/da.h"
#include "dret/backend/rlcsa.h"
#include "dret/backend/sa_phi.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_brute_backend.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_gcda_backend.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/doc_list/doc_list_slp.h"

#include "base_test.h"

namespace {

using S = std::reference_wrapper<dret::GenericStorage>;
using A = dret::Alphabet<>;
using RC = sri::RIndexCount<S, A>;

template <typename TIndex, typename... TArgs>
struct Kind {
  using Index = TIndex;
  static Index Make(S t_storage) { return Index(t_storage); }
};

struct BruteRI : Kind<dret::DocListIdxBrute<S, A, sri::RIndex<S>>> {};
struct BruteDA : Kind<dret::DocListIdxBruteBackend<S>> {};
struct BruteSAPhi : Kind<dret::DocListIdxBruteBackend<S, A, RC, dret::backend::SAPhi<S, A::int_width>>> {};
struct BruteRLCSA : Kind<dret::DocListIdxBruteBackend<S, A, RC, dret::backend::RLCSA<S, A::int_width>>> {};
struct BruteGCDA : Kind<dret::DocListIdxSLP<S>> {};
struct GCDA : Kind<dret::gcda::DocListIdxGCDA<S>> {
  static Index Make(S t_storage) { return Index(t_storage, 8, 4); }
};
struct GCDADiff : Kind<dret::dgcda::DocListIdxDGCDA<S>> {
  static Index Make(S t_storage) { return Index(t_storage, 8, 4); }
};
struct GCDAOverDA : Kind<dret::gcda::DocListIdxGCDABackend<S, A, RC, dret::backend::PlainDA<S, A::int_width>>> {
  static Index Make(S t_storage) { return Index(t_storage, 8, 4); }
};
struct PDL : Kind<dret::pdl::DocListIdxPDLPlain<S>> {
  static Index Make(S t_storage) { return Index(t_storage, 8, 4.0f); }
};
struct PDLOverRLCSA : Kind<dret::pdl::DocListIdxPDLPlain<S, A, RC, dret::backend::RLCSA<S, A::int_width>>> {
  static Index Make(S t_storage) { return Index(t_storage, 8, 4.0f); }
};
struct Sada : Kind<dret::rmq::DocListIdxRMQ<S, A, RC, dret::rmq::SadaCore<S>>> {};
struct SadaV : Kind<dret::rmq::DocListIdxRMQ<S, A, RC, dret::rmq::SadaValueCore<S>>> {};
struct Ilcp : Kind<dret::rmq::DocListIdxRMQ<S, A, RC, dret::rmq::IlcpCore<S>>> {};
struct IlcpV : Kind<dret::rmq::DocListIdxRMQ<S, A, RC, dret::rmq::IlcpValueCore<S>>> {};
struct Cilcp : Kind<dret::rmq::DocListIdxRMQ<S, A, RC, dret::rmq::CilcpCore<S>>> {};
struct CilcpV : Kind<dret::rmq::DocListIdxRMQ<S, A, RC, dret::rmq::CilcpValueCore<S>>> {};

using Kinds = testing::Types<BruteRI,
                             BruteDA,
                             BruteSAPhi,
                             BruteRLCSA,
                             BruteGCDA,
                             GCDA,
                             GCDADiff,
                             GCDAOverDA,
                             PDL,
                             PDLOverRLCSA,
                             Sada,
                             SadaV,
                             Ilcp,
                             IlcpV,
                             Cilcp,
                             CilcpV>;

template <typename TKind>
class StreamRoundTripTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    // Versions of one random reference, a few cut short: repetitive, with documents of different
    // lengths, so every family has runs, sampled nodes and lists to store.
    std::mt19937 rng(20261125);
    std::uniform_int_distribution<int> base(0, 3), pos(0, 199), coin(0, 9);
    std::string ref;
    for (int i = 0; i < 200; ++i) ref += "ACGT"[base(rng)];
    for (int d = 0; d < 20; ++d) {
      auto doc = ref;
      for (int k = 0; k < 5; ++k) doc[pos(rng)] = "ACGT"[base(rng)];
      if (coin(rng) < 3) doc = doc.substr(0, 100 + pos(rng) / 2);
      data_ += doc;
      data_ += '\1';
    }
    Init(data_);

    std::set<std::string> patterns;
    std::uniform_int_distribution<std::size_t> start(0, data_.size() - 1), len(1, 8);
    while (patterns.size() < 400) {
      auto p = data_.substr(start(rng), len(rng));
      if (p.find('\1') == std::string::npos) patterns.insert(p);
    }
    patterns_.assign(patterns.begin(), patterns.end());
  }

  template <typename TIndex>
  std::vector<std::vector<std::size_t>> Answers(const TIndex& t_index) const {
    std::vector<std::vector<std::size_t>> answers;
    for (const auto& p : patterns_) {
      std::vector<std::size_t> docs;
      t_index.Search(p, [&docs](std::size_t d) { docs.push_back(d); });
      std::sort(docs.begin(), docs.end());
      answers.push_back(std::move(docs));
    }
    return answers;
  }

  std::string data_;
  std::vector<std::string> patterns_;
};

TYPED_TEST_SUITE(StreamRoundTripTest, Kinds);

TYPED_TEST(StreamRoundTripTest, LoadsWhatItSerialized) {
  dret::GenericStorage built_storage;
  auto built = TypeParam::Make(std::ref(built_storage));
  construct(built, this->config_);

  std::stringstream stream;
  const auto written = sdsl::serialize(built, stream);
  ASSERT_EQ(written, stream.str().size());

  // The stream alone: a load that reopened a cache file would fail now.
  std::filesystem::remove_all(this->tmp_dir_);

  dret::GenericStorage loaded_storage;
  auto loaded = TypeParam::Make(std::ref(loaded_storage));
  loaded.load(stream);
  EXPECT_EQ(static_cast<std::size_t>(stream.tellg()), written) << "the load did not read the stream to its end";

  std::stringstream again;
  sdsl::serialize(loaded, again);
  EXPECT_TRUE(again.str() == stream.str()) << "serializing the loaded index wrote other bytes";

  EXPECT_EQ(this->Answers(loaded), this->Answers(built));
  EXPECT_EQ(sdsl::size_in_bytes(loaded), written);  // not built's: the RLCSA's serialize reads its cache files
  const auto sizes_built = built.GetSizeReport();
  const auto sizes_loaded = loaded.GetSizeReport();
  ASSERT_EQ(sizes_loaded.size(), sizes_built.size());
  for (std::size_t i = 0; i < sizes_built.size(); ++i) {
    EXPECT_EQ(sizes_loaded[i].name, sizes_built[i].name);
    EXPECT_EQ(sizes_loaded[i].bytes, sizes_built[i].bytes) << sizes_built[i].name;
  }
}

}  // namespace
