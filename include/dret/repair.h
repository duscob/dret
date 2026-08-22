//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/22/26.
//
// Choosing how to run irepair over a document array.
//
// The distribution ships the balanced build twice, under `bal/` and under
// `large/bal/`. They are the same algorithm; the second indexes the sequence
// with `long long` instead of `int` and takes an extra <MB> argument. Called
// kBal32 and kBal64 here, because "balanced" is what they have in common and the
// position width is what separates them.
//
// Given enough MB the two emit byte-identical .R/.C, so which one runs is an
// operational choice -- except past 2^31 elements, where only kBal64 is correct.
//
// <MB> is NOT a memory cap. It is the point at which irepair stops using a
// memory-lean pass and switches to the fast in-memory algorithm, and below that
// point it produces a DIFFERENT (valid, deterministic) grammar. So the value is
// part of what gets built, not merely how fast it is built.
//
// This header is deliberately free of the filesystem, subprocesses and the
// environment: the whole policy is PlanRePair, a pure function of the element
// count and the options, so it can be asserted directly. The measurements the
// policy encodes are in docs/bug_empty_da_page_big.md.
//

#pragma once

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>

namespace dret {

namespace repair {

// Which irepair binary to run.
//
// Both are RePair's *balanced* build -- upstream ships them as `bal/` and
// `large/bal/`, and the second is a copy of the first, not a different
// algorithm. They differ in one line of basics.h: the type used for sequence
// positions, `int` versus `typedef long long relong`. Hence the names: the
// number is the width of a position, NOT the compilation target (both compile
// -m64).
//
// Given enough <MB> the two emit byte-identical .R/.C, which is what makes
// choosing per collection safe -- see PlanRePair.
enum class Variant { kAuto, kBal32, kBal64 };

// Everything RunRePair needs beyond the file itself. Defaults auto-derive from
// the input, so `RunRePair(path)` remains correct for every collection; the
// fields exist so a driver can record and override the choice explicitly. This
// is deliberately plain data and reads no environment: <MB> changes the grammar,
// so the value used must be visible in the run's configuration rather than in
// the ambient process state.
struct Options {
  Variant variant = Variant::kAuto;
  // Unset: derive the value that guarantees the fast path for this input.
  std::optional<std::uintmax_t> mb;
  // Unset: no cap. Set below the derived value to trade the grammar for memory.
  std::optional<std::uintmax_t> mb_max;
};

// kBal32 (`bal/irepair`) walks the sequence with `int len` / `int i`, so it is correct
// only while the element count fits in one. INT_MAX is 2,147,483,647; stay under
// it with room to spare, because positions are formed as `len + k` in places and
// the alphabet grows by one symbol per rule.
inline constexpr std::uintmax_t kIntSafeElems = 2'000'000'000ull;

// kBal64 (`large/bal/irepair`) skips its memory-lean pass -- and so reproduces kBal32's
// grammar exactly -- as soon as
//
//     (len / 1024 / 1024) * 3 * sizeof(relong)  <=  MB
//
// with relong = long long. That is the literal condition in its prepare0() and
// its main loop, so the MB needed for the fast path is computable from the input
// size alone. One MB of slack covers the truncating division.
inline constexpr std::uintmax_t FastPathMB(std::uintmax_t t_elems) {
  return (t_elems / 1024u / 1024u) * 3u * sizeof(long long) + 1u;
}

// Peak RSS of the fast path, measured at ~28-31 bytes per element (kBal32 needs
// ~19). Diagnostic only: it explains a decision, it never makes one.
inline constexpr std::uintmax_t FastPathPeakMB(std::uintmax_t t_elems) {
  return (t_elems / 1024u / 1024u) * 31u + 1u;
}

// The resolved decision for one document array. Pure data, so the policy can be
// asserted without a filesystem, a subprocess, or an environment -- which is the
// point: which grammar you get is a function of these fields, and the mapping
// from (size, options) to fields was established by measurement, not by reading
// the irepair source alone.
struct Plan {
  bool use_64bit = false;
  // Only meaningful when use_64bit; the 32-bit build takes no second argument.
  std::optional<std::uintmax_t> mb;
  // True when mb sits below FastPathMB, i.e. irepair will take its memory-lean
  // pass and emit a DIFFERENT (valid, deterministic) grammar. Callers should say
  // so out loud rather than let a cache quietly become inconsistent.
  bool lean_path = false;
};

// Resolve the plan for an array of t_elems elements. Throws only for a request
// that cannot be honoured: forcing the 32-bit build onto an array it cannot
// index would silently wrap, which is the whole defect this selection prevents.
inline Plan PlanRePair(std::uintmax_t t_elems, const Options& t_options = {}) {
  const bool needs_64bit = t_elems > kIntSafeElems;

  Plan plan;
  plan.use_64bit = (t_options.variant == Variant::kBal64) ||
                   (t_options.variant == Variant::kAuto && needs_64bit);

  if (!plan.use_64bit && needs_64bit) {
    throw std::runtime_error(
        "RePair: this document array holds " + std::to_string(t_elems) + " elements, past the " +
        std::to_string(kIntSafeElems) +
        " that the 32-bit build indexes with `int`, but bal32 was requested explicitly. "
        "Use auto or bal64.");
  }

  if (!plan.use_64bit) {
    return plan;
  }

  const auto fast_path = FastPathMB(t_elems);
  auto mb = t_options.mb.value_or(fast_path);
  if (t_options.mb_max.has_value()) {
    mb = std::min(mb, *t_options.mb_max);
  }
  if (mb == 0) {
    // irepair rejects <MB> <= 0, and a cap of zero is a configuration mistake
    // rather than a request for minimal memory.
    throw std::runtime_error("RePair: the resolved <MB> is 0; give a positive value.");
  }

  plan.mb = mb;
  plan.lean_path = mb < fast_path;
  return plan;
}

}  // namespace repair

}  // namespace dret
