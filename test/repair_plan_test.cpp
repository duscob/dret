//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/22/26.
//
// Direct tests of dret::repair::PlanRePair -- the decision of which irepair
// binary to run over a document array, and with what <MB>.
//
// These assertions are not derived from reading irepair's source. They pin
// numbers that were established by measurement on 2026-08-22 and written up in
// docs/bug_empty_da_page_big.md:
//
//   * On a 5,997,971-element array the fast-path threshold is 120 MB. Running
//     at 120 or above reproduced kBal32's .R/.C byte for byte; 119, 100 and 64
//     each produced a different (valid) grammar; 32 and 8 did not finish in 90s.
//   * On the real `page` document array (115,163,007 elements) the derived
//     value 2617 reproduced kBal32's grammar, and capping to 1000 produced a
//     different one -- reproducibly the same different one across runs.
//
// The point of testing the planner in isolation is that the mapping from
// (size, options) to (binary, MB) decides *which grammar gets built*, so it has
// to survive refactoring without re-running a day of experiments.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstdint>
#include <optional>
#include <stdexcept>

#include "dret/repair.h"

namespace {

using dret::repair::FastPathMB;
using dret::repair::kIntSafeElems;
using dret::repair::Options;
using dret::repair::PlanRePair;
using dret::repair::Variant;

// Element counts of collections this project actually builds.
constexpr std::uintmax_t kPageElems = 115'163'007;    // real/small/page
constexpr std::uintmax_t kSyntheticElems = 5'997'971; // the measured sweep above
constexpr std::uintmax_t kEfaecGenElems = 3'913'901'169;

//~~~~~~~ variant selection

TEST(RePairPlan, SmallArrayUsesBal32AndTakesNoMBArgument) {
  const auto plan = PlanRePair(kPageElems);

  EXPECT_FALSE(plan.use_64bit);
  // The 32-bit build is a one-argument program; handing it an <MB> would make it
  // print usage and exit non-zero.
  EXPECT_FALSE(plan.mb.has_value());
  EXPECT_FALSE(plan.lean_path);
}

TEST(RePairPlan, ArrayPastTheIntLimitUsesBal64) {
  const auto plan = PlanRePair(kEfaecGenElems);

  EXPECT_TRUE(plan.use_64bit);
  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_FALSE(plan.lean_path);
}

TEST(RePairPlan, TheSwitchHappensExactlyAtTheIntSafeLimit) {
  EXPECT_FALSE(PlanRePair(kIntSafeElems).use_64bit);
  EXPECT_TRUE(PlanRePair(kIntSafeElems + 1).use_64bit);
}

// The limit is deliberately below INT_MAX: irepair forms positions as `len + k`
// in places, so the last representable element count is not a safe one.
TEST(RePairPlan, TheIntSafeLimitLeavesHeadroomBelowIntMax) {
  EXPECT_LT(kIntSafeElems, static_cast<std::uintmax_t>(2'147'483'647));
}

TEST(RePairPlan, Bal64CanBeForcedOnAnArrayThatDoesNotNeedIt) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal64});

  EXPECT_TRUE(plan.use_64bit);
  ASSERT_TRUE(plan.mb.has_value());
  // Forcing kBal64 must not silently change the grammar, so the
  // derived MB still has to clear the fast-path threshold.
  EXPECT_FALSE(plan.lean_path);
}

// The whole point of the selection. Wrapping `int` indices would yield a corrupt
// grammar that still looks like a grammar, so this refuses instead.
TEST(RePairPlan, ForcingBal32PastTheIntLimitThrows) {
  EXPECT_THROW(PlanRePair(kEfaecGenElems, {.variant = Variant::kBal32}), std::runtime_error);
}

TEST(RePairPlan, ForcingBal32BelowTheIntLimitIsFine) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal32});

  EXPECT_FALSE(plan.use_64bit);
  EXPECT_FALSE(plan.mb.has_value());
}

//~~~~~~~ derived <MB>

// Measured: 120 reproduced kBal32's grammar, 119 did not. The derived value must
// land on the safe side of that boundary.
TEST(RePairPlan, DerivedMBClearsTheMeasuredThresholdOnTheSweptArray) {
  const auto plan = PlanRePair(kSyntheticElems, {.variant = Variant::kBal64});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_GE(*plan.mb, 120u);
  EXPECT_EQ(*plan.mb, 121u);
  EXPECT_FALSE(plan.lean_path);
}

// Measured on the real page DA: the run at 2617 matched kBal32 byte for byte.
TEST(RePairPlan, DerivedMBMatchesTheValueVerifiedOnPage) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal64});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 2617u);
}

// ~87.5 GB. Recorded so that a change in the formula shows up as a diff here
// rather than as an OOM on the compute host days into a campaign.
TEST(RePairPlan, DerivedMBForEfaecGenIsTheDocumentedFigure) {
  const auto plan = PlanRePair(kEfaecGenElems);

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 89'569u);
}

TEST(RePairPlan, DerivedMBGrowsWithTheArray) {
  EXPECT_LT(FastPathMB(kPageElems), FastPathMB(kEfaecGenElems));
}

// An array smaller than one megabyte truncates to zero in irepair's own
// comparison; the plan must still pass a positive <MB>, since it rejects <= 0.
TEST(RePairPlan, TinyArrayStillGetsAPositiveMB) {
  const auto plan = PlanRePair(1000, {.variant = Variant::kBal64});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_GT(*plan.mb, 0u);
  EXPECT_FALSE(plan.lean_path);
}

//~~~~~~~ explicit overrides

TEST(RePairPlan, ExplicitMBIsHonoured) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb = 4096});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 4096u);
  EXPECT_FALSE(plan.lean_path);
}

// Measured: capping page to 1000 produced a different grammar. The flag is what
// lets the caller say so rather than let a cache go quietly inconsistent.
TEST(RePairPlan, ExplicitMBBelowTheThresholdIsFlaggedAsLean) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb = 1000});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 1000u);
  EXPECT_TRUE(plan.lean_path);
}

TEST(RePairPlan, CapBelowTheDerivedValueBitesAndIsFlagged) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb_max = 1000});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 1000u);
  EXPECT_TRUE(plan.lean_path);
}

TEST(RePairPlan, CapAboveTheDerivedValueDoesNothing) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb_max = 999'999});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 2617u);
  EXPECT_FALSE(plan.lean_path);
}

TEST(RePairPlan, CapAppliesToAnExplicitMBToo) {
  const auto plan =
      PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb = 4096, .mb_max = 2000});

  ASSERT_TRUE(plan.mb.has_value());
  EXPECT_EQ(*plan.mb, 2000u);
  EXPECT_TRUE(plan.lean_path);
}

// irepair exits with usage on <MB> <= 0, and a zero cap is a mistake in the
// caller's configuration rather than a request for minimal memory.
TEST(RePairPlan, ZeroMBIsRejected) {
  EXPECT_THROW(PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb = 0}), std::runtime_error);
  EXPECT_THROW(PlanRePair(kPageElems, {.variant = Variant::kBal64, .mb_max = 0}),
               std::runtime_error);
}

// <MB> is meaningless to kBal32, so options that only concern it must
// not leak into a kBal32 invocation.
TEST(RePairPlan, MBOptionsAreIgnoredWhenBal32IsChosen) {
  const auto plan = PlanRePair(kPageElems, {.variant = Variant::kBal32, .mb = 4096, .mb_max = 8});

  EXPECT_FALSE(plan.use_64bit);
  EXPECT_FALSE(plan.mb.has_value());
  EXPECT_FALSE(plan.lean_path);
}

//~~~~~~~ defaults

// `RunRePair(path)` has to stay correct for every collection without a caller
// supplying anything, which is what makes the options purely additive.
TEST(RePairPlan, DefaultOptionsNeverSelectTheLeanPath) {
  for (const auto elems : {std::uintmax_t{1}, std::uintmax_t{1'000'000}, kSyntheticElems,
                           kPageElems, kIntSafeElems, kEfaecGenElems}) {
    const auto plan = PlanRePair(elems);
    EXPECT_FALSE(plan.lean_path) << "elems = " << elems;
    if (plan.use_64bit) {
      ASSERT_TRUE(plan.mb.has_value()) << "elems = " << elems;
      EXPECT_GE(*plan.mb, FastPathMB(elems)) << "elems = " << elems;
    }
  }
}

}  // namespace
