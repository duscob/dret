//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/13/26.
//
// DA / GCDA / DGCDA raw-range policy parity tests, Task 38 of
// docs/pdl_indexes_tasks.md. The three PDLRawRangePolicy backings —
// PDLGetDocsDA (plain int_vector), PDLGetDocsSLP (LightSLP-compressed),
// PDLGetDocsDSLP (DifferentialLightSLP-compressed) — all expose the
// same getDocs(sp, ep, report) contract. This test asserts they
// produce identical document-id sequences for every probed [sp, ep)
// range against a small real dataset.
//
// Acceptance: DA, GCDA, and DGCDA policies report identical document
// ids for every tested range.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "dret/config.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"
#include "dret/rmq/rmq_get_doc_policies.h"

#include "base_test.h"

namespace {

using ExternalGenericStorage = std::reference_wrapper<dret::GenericStorage>;

class PDLGetDocsPolicyParityTest : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
    // Bootstrap a PDL Plain index whose construct() builds Text / SA /
    // DocEnds / DA / LCP caches. The standalone get-doc policies below
    // can then construct() against the existing DA cache (DA is no-op,
    // SLP/DSLP construct their LightSLP/DifferentialLightSLP entries).
    dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage> bootstrap(
        std::ref(storage_));
    construct(bootstrap, config_);
  }

  sri::GenericStorage storage_;
  // 5-doc fixture (same as pdl_search_test.cpp): exercises each get-doc
  // policy over a non-trivial DA without making the test slow.
  const std::string data_ =
      "abracadabra\1barbara\1albatross\1amanda\1lasagna\1";
};

// Collect getDocs(sp, ep) into a sorted vector. The PDL search path
// collects doc ids and sort-dedupes before reporting, so the contract
// the three policies must share is multi-set equality, not order
// equality. (Internally the SLP expansion descends right-to-left when
// the start position does not align to a leaf boundary, so the raw
// stream order differs from DA's straight forward iteration. That's
// invisible to PDL callers but breaks an order-sensitive comparison.)
template <typename TPolicy>
std::vector<std::size_t> CollectGetDocsSorted(const TPolicy& policy,
                                              std::size_t sp,
                                              std::size_t ep) {
  std::vector<std::size_t> out;
  auto report = [&out](std::size_t d) { out.push_back(d); };
  policy.getDocs(sp, ep, report);
  std::sort(out.begin(), out.end());
  return out;
}

TEST_F(PDLGetDocsPolicyParityTest, DA_SLP_DSLP_ReportIdenticalDocsAcrossRanges) {
  using DAPolicy   = dret::pdl::PDLGetDocsDA<ExternalGenericStorage>;
  using SLPPolicy  = dret::pdl::PDLGetDocsSLP<ExternalGenericStorage>;
  using DSLPPolicy = dret::pdl::PDLGetDocsDSLP<ExternalGenericStorage>;

  // Build each policy with the same external storage. The PDL adapter
  // has no ctor accepting storage directly, so wrap an explicitly-
  // constructed inner. SLP/DSLP take (storage, block_size,
  // storing_factor); their defaults (512, 4) match the bootstrap's.
  DAPolicy   da{   typename DAPolicy::Inner  {std::ref(storage_)} };
  SLPPolicy  slp{  typename SLPPolicy::Inner {std::ref(storage_), 512, 4.0f} };
  DSLPPolicy dslp{ typename DSLPPolicy::Inner{std::ref(storage_), 512, 4.0f} };

  // Construct: DA is a no-op (DA cache already exists). SLP/DSLP build
  // their LightSLP / DifferentialLightSLP cache entries from the DA.
  construct(da, config_);
  construct(slp, config_);
  construct(dslp, config_);

  // Load each policy's inner so it sees the cached structure.
  da.inner().load(config_);
  slp.inner().load(config_);
  dslp.inner().load(config_);

  const std::size_t n = da.inner().size();
  ASSERT_GT(n, 0u);

  // Probe ranges: empty, full, prefix, suffix, mid-slice, single
  // positions across the SA. Cover every shape that getDocs is
  // expected to handle.
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 0},               // empty
      {0, n},               // full
      {0, 1},               // single first
      {n - 1, n},           // single last
      {0, n / 2},           // prefix
      {n / 2, n},           // suffix
      {n / 4, 3 * n / 4},   // mid-slice
      {1, n - 1},           // strict interior
      {2, 7},               // arbitrary small
      {5, 5},               // empty mid
  };

  for (auto [sp, ep] : ranges) {
    SCOPED_TRACE(::testing::Message() << "[" << sp << ", " << ep << ")");
    auto da_out   = CollectGetDocsSorted(da, sp, ep);
    auto slp_out  = CollectGetDocsSorted(slp, sp, ep);
    auto dslp_out = CollectGetDocsSorted(dslp, sp, ep);

    // Element count must also match — catches a policy emitting too
    // many or too few docs even if multiset equality happened to
    // coincide on a small fixture.
    ASSERT_EQ(da_out.size(), ep - sp);
    EXPECT_EQ(slp_out.size(), da_out.size());
    EXPECT_EQ(dslp_out.size(), da_out.size());

    EXPECT_EQ(slp_out, da_out)  << "GCDA-SLP disagrees with DA";
    EXPECT_EQ(dslp_out, da_out) << "DGCDA disagrees with DA";
  }
}

}  // namespace
