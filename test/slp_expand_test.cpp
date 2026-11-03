//
// ExpandSLP over every GCDA representation against the document array: ranges inside one sampled
// leaf (the nearer end on either side), ranges ending or starting on a leaf boundary, ranges over
// several leaves. The skip that reaches a mid-leaf position is the part these representations do
// differently (slp_tools.h, compact_louds_slp.h), so it is what this test is for.
//

#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "dret/doc_list/doc_list_gcda.h"
#include "dret/slp/slp_tools.h"

#include "base_test.h"

namespace {

using S = std::reference_wrapper<dret::GenericStorage>;
using RC = sri::RIndexCount<S, dret::Alphabet<>>;
using SLPCNF = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
using Light = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>, grammar::SampledSLP<>,
                                grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
using Combined = grammar::CombinedSLPWithUnitCover<grammar::CombinedSLP<SLPCNF, grammar::SampledSLP<>, sdsl::int_vector<>>>;
using SetsPlain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;

// Which representations skip by stored span lengths (slp_tools.h SkipWhole).
static_assert(dret::StoresSpanLength<Combined, std::size_t>);
static_assert(!dret::StoresSpanLength<Light, std::size_t>);
static_assert(!dret::StoresSpanLength<grammar::CompactBPSLP<>, std::size_t>);
static_assert(!dret::StoresSpanLength<grammar::CompactLOUDSSLP<>, std::size_t>);

// The index's grammar, which DocListIdxGCDA keeps protected.
template <typename TSLP>
struct Probe : dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, TSLP, SetsPlain> {
  using Base = dret::gcda::DocListIdxGCDA<S, dret::Alphabet<>, RC, TSLP, SetsPlain>;
  using Base::Base;
  const TSLP& slp() const { return *this->slp_; }
};

class SLPExpandTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    std::mt19937 rng(20261103);
    std::uniform_int_distribution<int> base(0, 3), pos(0, 239), coin(0, 9);
    std::string ref;
    for (int i = 0; i < 240; ++i) ref += "ACGT"[base(rng)];
    std::string data;
    for (int d = 0; d < 24; ++d) {
      auto doc = ref;
      for (int k = 0; k < 6; ++k) doc[pos(rng)] = "ACGT"[base(rng)];
      if (coin(rng) < 3) doc = doc.substr(0, 120 + pos(rng) / 2);
      data += doc;
      data += '\1';
    }
    Init(data);
  }

  // Every range kind at every start, plus random ranges; the expansion must be DA[bp, ep) (order free:
  // the backward paths report it reversed, and document listing deduplicates).
  template <typename TSLP>
  void Check(std::uint32_t t_block_size) {
    Probe<TSLP> index(std::ref(storage_), t_block_size, 4);
    construct(index, config_);
    sdsl::int_vector<> da;
    ASSERT_TRUE(sdsl::load_from_cache(da, config_.keys[dret::conf::kDA].get<std::string>(), config_, true));
    const std::size_t n = da.size();
    const auto& slp = index.slp();

    std::vector<std::pair<std::size_t, std::size_t>> ranges;
    const std::size_t b = t_block_size;
    for (std::size_t bp = 0; bp < n; ++bp)
      for (std::size_t len : {std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{7}, b / 2, b - 1, b, b + 1, 2 * b + 3})
        if (len && bp + len <= n) ranges.emplace_back(bp, bp + len);
    std::mt19937 rng(t_block_size);
    std::uniform_int_distribution<std::size_t> any(0, n);
    for (int i = 0; i < 4000; ++i) {
      auto x = any(rng), y = any(rng);
      if (x != y) ranges.emplace_back(std::min(x, y), std::max(x, y));
    }

    std::size_t wrong = 0;
    for (const auto& [bp, ep] : ranges) {
      std::vector<std::size_t> got;
      auto report = [&got](auto v) { got.emplace_back(static_cast<std::size_t>(v)); };
      dret::ExpandSLP(slp, bp, ep, report);
      std::vector<std::size_t> want(da.begin() + bp, da.begin() + ep);
      std::sort(got.begin(), got.end());
      std::sort(want.begin(), want.end());
      if (got != want && wrong++ < 5)
        ADD_FAILURE() << "b=" << b << " [" << bp << ", " << ep << "): " << got.size() << " values, want " << want.size();
    }
    EXPECT_EQ(wrong, 0u) << "b=" << b << ", " << ranges.size() << " ranges";
  }

  dret::GenericStorage storage_;
};

TEST_F(SLPExpandTest, cached) {
  for (std::uint32_t b : {4u, 16u, 64u, 256u}) Check<Light>(b);
}

TEST_F(SLPExpandTest, ondemand) {
  for (std::uint32_t b : {4u, 16u, 64u, 256u}) Check<Combined>(b);
}

TEST_F(SLPExpandTest, compact_bp) {
  for (std::uint32_t b : {4u, 16u, 64u, 256u}) Check<grammar::CompactBPSLP<>>(b);
}

TEST_F(SLPExpandTest, compact_louds) {
  for (std::uint32_t b : {4u, 16u, 64u, 256u}) Check<grammar::CompactLOUDSSLP<>>(b);
}

}  // namespace
