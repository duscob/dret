//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// The grammars of the DA without a sampled tree, as Brute-GCDA (doc_list/doc_list_slp.h) and the
// bare-grammar backend (backend/grammar.h) read them: the CNF grammar alone, the base differential
// grammar sampled every few positions, or a sampled GCDA grammar without its document lists. For
// each, SlpNsTraits names its cache components and builds them, and SlpNsExpandRange expands a range
// of it.
//

#pragma once

#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/memory_management.hpp>

#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "dret/cache_components.h"
#include "dret/config.h"
#include "dret/slp/da_grammar.h"
#include "dret/slp/differential_slp.h"
#include "dret/slp/gcda_grammar.h"
#include "dret/slp/slp_components.h"
#include "dret/slp/slp_tools.h"

namespace dret {

// Forward decl of the plain bare-SLP build (defined below).
template <typename TVarsContainer, typename TLengthsContainer>
void construct(grammar::SLP<TVarsContainer, TLengthsContainer>& t_slp,
               Config& t_config, const std::string& t_datafile);

// Default sample spacing of the non-sampled differential SLP (bare-diff), the
// one GCDA-nolists has always used. Differential decode needs *some* anchors;
// the spacing is the only knob of that index.
inline constexpr std::uint32_t kDiffBlockSize = 512;

// Storing factor of a sampled GCDA grammar without its document lists. Besides
// the lists, which this index does not load, it decides which nodes enter the
// sampled tree (AreChildrenTooBig), so it also shapes the tree this index
// expands through, slightly (page, compact-louds at block size 1024: 2.77 MB at
// 4, 2.73 MB at inf). 32 is the one the RMQ/PDL GCDA backend uses.
inline constexpr float kNolistsStoringFactor = 32.0f;

//~~~~~~~  SLP-NS range expansion, dispatched by SLP type  ~~~~~~~

// Plain grammar::SLP: split [sp, ep) into a span cover, expand each variable.
template <typename TVarsContainer, typename TLengthsContainer, typename TReport>
void SlpNsExpandRange(const grammar::SLP<TVarsContainer, TLengthsContainer>& slp,
                      std::size_t sp, std::size_t ep, TReport add) {
  std::vector<typename grammar::SLP<TVarsContainer, TLengthsContainer>::VariableType> cover;
  grammar::ComputeSpanCover(slp, sp, ep, std::back_inserter(cover));
  for (auto var : cover) {
    auto length = slp.SpanLength(var);
    grammar::ExpandSLPForward(slp.GetRules(), slp.Sigma(), var, length, add);
  }
}

// Base differential SLP: range-decompress [sp, ep) directly (no sampled tree).
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples,
          typename TSampleRootsPos, typename TBV, typename TReport>
void SlpNsExpandRange(const DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& slp,
                      std::size_t sp, std::size_t ep, TReport add) {
  ExpandSLP(slp, sp, ep, add);
}

// Sampled GCDA grammar: expand [sp, ep) through its sampled tree, as the RMQ/PDL
// GCDA backend does. The two overloads above are more specialised and win.
template <typename TSLP, typename TReport>
void SlpNsExpandRange(const TSLP& slp, std::size_t sp, std::size_t ep, TReport add) {
  ExpandSLP(slp, sp, ep, add);
}

//~~~~~~~  SLP-NS cache key + build, dispatched by SLP type  ~~~~~~~

// Sampled GCDA grammar: the components of the GCDA cell (block size = the
// spacing, storing factor) that the RMQ/PDL GCDA backend loads, shared with it.
template <typename TSLP>
struct SlpNsTraits {
  static std::vector<Component> components(const JSON& t_keys, uint32_t t_spacing,
                                           float t_storing_factor = kNolistsStoringFactor) {
    return CacheComponents(TSLP{}, t_keys, SampledTreeCell{t_spacing, t_storing_factor});
  }
  static void build(Config& t_config, uint32_t t_spacing, float t_storing_factor = kNolistsStoringFactor) {
    const auto comps = components(t_config.keys, t_spacing, t_storing_factor);
    if (ComponentsExist(comps, t_config)) return;
    auto event = sdsl::memory_monitor::event(comps.back().key);
    auto da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[conf::kDA].get<std::string>(), t_config);
    TSLP slp;
    gcda::construct(slp, t_config, da, t_spacing, t_storing_factor);
  }
};

// Plain bare grammar::SLP: the CNF grammar of the DA alone, the component the
// on-demand GCDA representation also stores. It has no sampling, so it ignores
// the spacing.
template <typename TVarsContainer, typename TLengthsContainer>
struct SlpNsTraits<grammar::SLP<TVarsContainer, TLengthsContainer>> {
  using TSLP = grammar::SLP<TVarsContainer, TLengthsContainer>;
  static std::vector<Component> components(const JSON& t_keys, uint32_t /*t_spacing*/, float /*t_storing_factor*/ = 0) {
    return CacheComponents(TSLP{}, t_keys);
  }
  static void build(Config& t_config, uint32_t t_spacing, float /*t_storing_factor*/ = 0) {
    if (ComponentsExist(components(t_config.keys, t_spacing), t_config)) return;
    auto da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[conf::kDA].get<std::string>(), t_config);
    TSLP slp;
    construct(slp, t_config, da);
  }
};

// Bare-diff — base DifferentialSLP sampled every t_spacing positions; every
// spacing shares the grammar, roots and span sums, and stores its own samples.
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples,
          typename TSampleRootsPos, typename TBV>
struct SlpNsTraits<DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>> {
  using TDiff = DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>;
  static std::vector<Component> components(const JSON& t_keys, uint32_t t_spacing, float /*t_storing_factor*/ = 0) {
    return CacheComponents(TDiff{}, t_keys, t_spacing);
  }
  static void build(Config& t_config, uint32_t t_spacing, float /*t_storing_factor*/ = 0) {
    if (ComponentsExist(components(t_config.keys, t_spacing), t_config)) return;
    auto event = sdsl::memory_monitor::event(components(t_config.keys, t_spacing).back().key);
    TDiff dslp;
    construct(dslp, t_config, t_spacing);
  }
};

//~~~~~~~


template <typename TVarsContainer, typename TLengthsContainer>
void construct(grammar::SLP<TVarsContainer, TLengthsContainer>& t_slp,
               Config& t_config,
               const std::string& t_datafile) {
  using namespace conf;

  // Build into a default-typed SLP first (push_back-friendly std::vector<uint32_t>
  // containers), then convert via the SLP template copy constructor — sdsl::dac_vector
  // and sdsl::vlc_vector have no push_back, so RePairReader cannot fill them
  // directly. The conversion uses grammar::Construct overloads (grammar/utility.h)
  // which dispatch via SFINAE: dac_vector / vlc_vector use their templated Container
  // constructor (which bit-packs internally); sdsl::int_vector falls through to the
  // resize-and-std::copy branch (which leaves the default 64-bit width). Pass a
  // bit-compress action so int_vector targets get tight per-element width; the
  // action no-ops on dac_vector / vlc_vector via `if constexpr` detection.
  const auto tmp_slp = LoadOrBuildDaCnfGrammar(t_config, t_datafile);
  auto bc = []<typename TVec>(TVec& v) {
    // Probe for bit_resize/width members directly — `requires { bit_compress(vv); }`
    // would unhelpfully pass for any T (since bit_compress is an unconstrained
    // template) and only fail at body instantiation. int_vector has both members;
    // std::vector and dac_vector / vlc_vector have neither.
    if constexpr (requires(TVec& vv) {
                    vv.bit_resize(std::size_t{});
                    vv.width(uint8_t{});
                  }) {
      sdsl::util::bit_compress(v);
    }
  };
  t_slp = grammar::SLP<TVarsContainer, TLengthsContainer>(tmp_slp, bc, bc);

  StoreComponents(t_slp, CacheComponents(t_slp, t_config.keys), t_config);
}

}  // namespace dret
