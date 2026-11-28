//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Document-array backends over the grammars of GCDA (see backend/backend.h for the contract):
//
//   GCDAGrammar      the sampled GCDA grammar (grammar::LightSLP and the other sampled representations),
//                    the same components GCDA loads
//   BareGrammar      a GCDA grammar without its sampled tree: the CNF grammar alone, or the differential
//                    grammar; the backend of Brute-GCDA
//   DiffGCDAGrammar  the sampled differential grammar (DifferentialLightSLP), as GCDA-differential loads it
//
// Each loads the components the corresponding index caches, so an index over one of them and that
// index share the files.

#pragma once

#include <cstddef>
#include <format>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "dret/backend/backend.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"
#include "dret/slp/bare_grammar.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/slp/differential_slp.h"
#include "dret/slp/gcda_grammar.h"
#include "dret/slp/slp_tools.h"

namespace dret {

// Trait: is T an instantiation of dret::DifferentialSLP<...>?  Used by
// BareGrammar::construct to switch construct overloads when the chosen
// bare-SLP type is a DifferentialSLP (bare-diff) rather than a plain
// grammar::SLP.
template <typename T>
struct is_differential_slp : std::false_type {};
template <typename TSLP, typename TRoots, typename TSpanSums,
          typename TSamples, typename TSampleRootsPos, typename TBV>
struct is_differential_slp<DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>>
    : std::true_type {};
template <typename T>
inline constexpr bool is_differential_slp_v = is_differential_slp<T>::value;

// ExpandSLP overload for bare `grammar::SLP<>` (no sampled tree, no covers).
// The default templated `ExpandSLP` (slp_tools.h) calls `Leaf`/`Position`/`Cover`,
// which `grammar::SLP<>` does not provide. This overload instead splits [bp, ep)
// into cover variables via `ComputeSpanCover` and expands each via
// `ExpandSLPForward` — same algorithm as `DocListIdxSLP::Search`. Must be visible
// before `BareGrammar`'s definition for ordinary two-phase lookup to find it.
template <typename TVars, typename TLens, typename Report>
void ExpandSLP(const grammar::SLP<TVars, TLens>& slp,
               std::size_t bp, std::size_t ep, Report& report) {
  if (bp >= ep) return;
  using V = typename grammar::SLP<TVars, TLens>::VariableType;
  std::vector<V> cover;
  grammar::ComputeSpanCover(slp, bp, ep, std::back_inserter(cover));
  for (auto var : cover) {
    auto length = slp.SpanLength(var);
    grammar::ExpandSLPForward(slp.GetRules(), slp.Sigma(), var, length, report);
  }
}

// In-order (the cover is left to right) but without a way to stop inside
// grammar::ExpandSLPForward, so this drops what follows a stop request.
template <typename TVars, typename TLens, typename Report>
void ExpandSLPUntil(const grammar::SLP<TVars, TLens>& slp, std::size_t bp, std::size_t ep, Report& report) {
  bool going = true;
  auto filtered = [&report, &going](auto d) {
    if (going)
      going = report(d);
  };
  ExpandSLP(slp, bp, ep, filtered);
}

}  // namespace dret

namespace dret::backend {

// Document-array lookup over a grammar-compressed LightSLP. By default the SLP
// type exactly matches gcda::DocListIdxGCDA's default TSLP so both indexes share
// the same typed cache entry.
template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>>
class GCDAGrammar : public FieldsComponent<GCDAGrammar<TStorage, t_width, TSLP>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<GCDAGrammar, TStorage, t_width>;
  using SLP = TSLP;
  using typename Base::size_type;

  explicit GCDAGrammar(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  GCDAGrammar() = default;

  std::size_t operator()(std::size_t i) const {
    std::size_t value = 0;
    auto report = [&value](auto d) {
      value = static_cast<std::size_t>(d);
    };
    ExpandSLP(*slp_, i, i + 1, report);
    return value;
  }

  // Half-open range expansion: calls r(doc) for each position in [b, e).
  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& r) const {
    ExpandSLP(*slp_, b, e, r);
  }

  // In-order expansion of [b, e) that stops once f returns false.
  template <typename F>
  void ExpandUntil(std::size_t b, std::size_t e, F&& f) const {
    auto report = [&f](auto d) { return f(static_cast<std::size_t>(d)); };
    ExpandSLPUntil(*slp_, b, e, report);
  }

  uint32_t block_size() const {
    return block_size_;
  }

  float storing_factor() const {
    return storing_factor_;
  }

  // The grammar, with its sampled tree: GCDA covers a range on it (sampled_tree_core.h, GrammarTree).
  const TSLP& slp() const { return *slp_; }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto components = [&t_self](const JSON& k) {
      return CacheComponents(TSLP{}, k, SampledTreeCell{t_self.block_size_, t_self.storing_factor_});
    };
    return std::tuple{Components<SizeAs::kCollect>("slp", t_self.slp_, components, "slp_")};
  }

  const TSLP* slp_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

template <typename TStorage, uint8_t t_width, typename TSLP>
void construct(GCDAGrammar<TStorage, t_width, TSLP>& t_get_doc, Config& t_config) {
  using namespace conf;

  const auto components =
      CacheComponents(TSLP{}, t_config.keys, SampledTreeCell{t_get_doc.block_size(), t_get_doc.storing_factor()});
  if (ComponentsExist(components, t_config))
    return;

  auto event = sdsl::memory_monitor::event(components.back().key);
  auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
  TSLP slp;
  // Every sampled representation, the differential one included (it reads the DA itself).
  gcda::GCDAVariantTraits<TSLP>::constructSLP(
      slp, t_config, filepath_da, t_get_doc.block_size(), t_get_doc.storing_factor());
}

// Document-array lookup over a GCDA grammar without its document lists: the CNF
// grammar alone (a plain grammar::SLP), the differential grammar without the
// sampled tree, or a sampled GCDA grammar whose lists are not loaded. The backend
// of Brute-GCDA (DocListIdxSLP). The block size is the sample spacing of a
// differential TSLP (a lookup jumps to a sample before skipping by span lengths)
// or the block size of a sampled one, with the storing factor; a plain TSLP has
// no sampling and ignores both (SlpNsTraits).
template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>
class BareGrammar : public FieldsComponent<BareGrammar<TStorage, t_width, TSLP>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<BareGrammar, TStorage, t_width>;
  using SLP = TSLP;
  using typename Base::size_type;

  explicit BareGrammar(const TStorage& t_storage,
                       uint32_t t_block_size = kDiffBlockSize,
                       float t_storing_factor = kNolistsStoringFactor)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  BareGrammar() = default;

  uint32_t block_size() const {
    return block_size_;
  }

  float storing_factor() const { return storing_factor_; }

  std::size_t operator()(std::size_t i) const {
    std::size_t value = 0;
    auto report = [&value](auto d) {
      value = static_cast<std::size_t>(d);
    };
    ExpandSLP(*slp_, i, i + 1, report);
    return value;
  }

  // Half-open range expansion: calls r(doc) for each position in [b, e).
  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& r) const {
    ExpandSLP(*slp_, b, e, r);
  }

  // In-order expansion of [b, e) that stops once f returns false.
  template <typename F>
  void ExpandUntil(std::size_t b, std::size_t e, F&& f) const {
    auto report = [&f](auto d) { return f(static_cast<std::size_t>(d)); };
    ExpandSLPUntil(*slp_, b, e, report);
  }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto components = [&t_self](const JSON& k) {
      return SlpNsTraits<TSLP>::components(k, t_self.block_size_, t_self.storing_factor_);
    };
    return std::tuple{Components("slp", t_self.slp_, components)};
  }

  const TSLP* slp_ = nullptr;
  uint32_t block_size_ = kDiffBlockSize;
  float storing_factor_ = kNolistsStoringFactor;
};

// Build the grammar's components. The DA they are built from is the index's: its
// construct() ensures it (EnsureDA, dret/collection.h) before it constructs the backend.
template <typename TStorage, uint8_t t_width, typename TSLP>
void construct(BareGrammar<TStorage, t_width, TSLP>& t_get_doc, Config& t_config) {
  SlpNsTraits<TSLP>::build(t_config, t_get_doc.block_size(), t_get_doc.storing_factor());
}

// Document-array lookup over a differential grammar-compressed SLP. The default
// type matches dgcda::DocListIdxDGCDA's default TSLP for typed cache sharing.
template <typename TStorage = GenericStorage, uint8_t t_width = 8, typename TDSLP = DifferentialLightSLP<>>
class DiffGCDAGrammar : public FieldsComponent<DiffGCDAGrammar<TStorage, t_width, TDSLP>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<DiffGCDAGrammar, TStorage, t_width>;
  using DSLP = TDSLP;
  using typename Base::size_type;

  explicit DiffGCDAGrammar(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  DiffGCDAGrammar() = default;

  std::size_t operator()(std::size_t i) const {
    std::size_t value = 0;
    auto report = [&value](auto d) {
      value = static_cast<std::size_t>(d);
    };
    ExpandSLP(*dslp_, i, i + 1, report);
    return value;
  }

  // Half-open range expansion: calls r(doc) for each position in [b, e).
  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& r) const {
    ExpandSLP(*dslp_, b, e, r);
  }

  // In-order expansion of [b, e) that stops once f returns false.
  template <typename F>
  void ExpandUntil(std::size_t b, std::size_t e, F&& f) const {
    auto report = [&f](auto d) { return f(static_cast<std::size_t>(d)); };
    ExpandSLPUntil(*dslp_, b, e, report);
  }

  uint32_t block_size() const {
    return block_size_;
  }

  float storing_factor() const {
    return storing_factor_;
  }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto components = [&t_self](const JSON& k) {
      return CacheComponents(TDSLP{}, k, SampledTreeCell{t_self.block_size_, t_self.storing_factor_});
    };
    return std::tuple{Components<SizeAs::kCollect>("dslp", t_self.dslp_, components, "dslp_")};
  }

  const TDSLP* dslp_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

template <typename TStorage, uint8_t t_width, typename TDSLP>
void construct(DiffGCDAGrammar<TStorage, t_width, TDSLP>& t_get_doc, Config& t_config) {
  using namespace conf;

  const auto components =
      CacheComponents(TDSLP{}, t_config.keys, SampledTreeCell{t_get_doc.block_size(), t_get_doc.storing_factor()});
  if (ComponentsExist(components, t_config))
    return;

  auto event = sdsl::memory_monitor::event(components.back().key);
  TDSLP dslp;
  dret::construct(dslp, t_config, t_get_doc.block_size(), t_get_doc.storing_factor());
}

}  // namespace dret::backend
