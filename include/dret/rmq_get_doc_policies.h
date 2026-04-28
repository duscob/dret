//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// GetDoc policies for RMQ document-listing indexes. Each policy implements
// single-position and range-based document-array lookup via a common interface:
//
//   std::size_t operator()(std::size_t i) const     — doc at SA position i
//   void operator()(std::size_t b, std::size_t e, TReport&) const  — half-open [b,e)
//   std::size_t size() const                        — total SA length, when cheaply available
//
// GetDocDA: plain sdsl::int_vector<> document array (no grammar compression).
// GetDocSLP: LightSLP document array shared with GCDA.
// GetDocDSLP: DifferentialLightSLP document array shared with DGCDA.

#pragma once

#include <cstddef>
#include <format>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "config.h"
#include "doc_list_idx_slp.h"
#include "doc_list_sampled_tree_dgcda.h"
#include "doc_list_sampled_tree_gcda.h"
#include "index_base.h"
#include "size_report.h"
#include "slp_tools.h"

namespace dret {

// ExpandSLP overload for bare `grammar::SLP<>` (no sampled tree, no covers).
// The default templated `ExpandSLP` (slp_tools.h) calls `Leaf`/`Position`/`Cover`,
// which `grammar::SLP<>` does not provide. This overload instead splits [bp, ep)
// into cover variables via `ComputeSpanCover` and expands each via
// `ExpandSLPForward` — same algorithm as `DocListIdxSLP::Search`. Must be visible
// before `GetDocSLP_NS`'s definition for ordinary two-phase lookup to find it.
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

}  // namespace dret

namespace dret::rmq {

// Document-array lookup over a plain sdsl::int_vector<>.
// Loads from the cache entry at conf::kDA, shared with all RMQ core constructors.
template <typename TStorage = GenericStorage, uint8_t t_width = 8>
class GetDocDA : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit GetDocDA(const TStorage& t_storage) : Base(t_storage) {}

  GetDocDA() = default;

  std::size_t operator()(std::size_t i) const {
    return static_cast<std::size_t>((*da_)[i]);
  }

  // Half-open range expansion: calls r(doc) for each position in [b, e).
  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& r) const {
    for (std::size_t i = b; i < e; ++i)
      r(static_cast<std::size_t>((*da_)[i]));
  }

  // Total number of SA positions. Used by IlcpLikeCore for the last-run sentinel.
  std::size_t size() const {
    return da_ ? da_->size() : 0;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += da_ ? sdsl::serialize(*da_, out, child, "da")
                   : sdsl::serialize_empty_object<sdsl::int_vector<>>(out, child, "da");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (da_)
      append(r, "da", sdsl::size_in_bytes(*da_));
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    key_da_ = t_keys[conf::kDA].get<std::string>();
    da_ = this->template loadItemPtr<sdsl::int_vector<>>(key_da_, t_source, true);
  }

  std::string key_da_;
  const sdsl::int_vector<>* da_ = nullptr;
};

// No-op: the DA cache entry is built by EnsureBasicStructures inside the core's
// construct() call before this would ever be invoked.
template <typename TStorage, uint8_t t_width>
void construct(GetDocDA<TStorage, t_width>& /*unused*/, Config& /*unused*/) {}

// Document-array lookup over a grammar-compressed LightSLP. By default the SLP
// type exactly matches gcda::DocListIdxGCDA's default TSLP so both indexes share
// the same typed cache entry.
template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>>
class GetDocSLP : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using SLP = TSLP;
  using typename Base::size_type;

  explicit GetDocSLP(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  GetDocSLP() = default;

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

  uint32_t block_size() const {
    return block_size_;
  }

  float storing_factor() const {
    return storing_factor_;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return slp_ ? sdsl::serialize(*slp_, out, child, "slp") : sdsl::serialize_empty_object<TSLP>(out, child, "slp");
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (slp_)
      collectSizes(r, *slp_, "slp_");
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    key_slp_ = std::format("{}-{}_", block_size_, storing_factor_) + t_keys[conf::kGCDA][conf::kSLP].get<std::string>();
    slp_ = this->template loadItemPtr<TSLP>(key_slp_, t_source, true);
  }

  std::string key_slp_;
  const TSLP* slp_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

template <typename TStorage, uint8_t t_width, typename TSLP>
void construct(GetDocSLP<TStorage, t_width, TSLP>& t_get_doc, Config& t_config) {
  using namespace conf;

  const auto key_prefix = std::format("{}-{}_", t_get_doc.block_size(), t_get_doc.storing_factor());
  const auto key_slp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
  if (sdsl::cache_file_exists<TSLP>(key_slp, t_config))
    return;

  auto event = sdsl::memory_monitor::event(key_slp);
  auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
  TSLP slp;
  gcda::construct(slp, t_config, filepath_da, t_get_doc.block_size(), t_get_doc.storing_factor());
}

// Document-array lookup over a bare `grammar::SLP<>` (Phase C cache, kSLPNS).
// Distinct from `GetDocSLP`: the bare SLP has no sampled tree, no precomputed
// covers, and no block_size / storing_factor knobs — its cache file is keyed
// directly on `kSLPNS` with no `{bs}-{sf}_` prefix, and is shared with
// `dret::DocListIdxSLP` via SDSL type-hashing on `grammar::SLP<>`.
template <typename TStorage = GenericStorage, uint8_t t_width = 8, typename TSLP = grammar::SLP<>>
class GetDocSLP_NS : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using SLP = TSLP;
  using typename Base::size_type;

  explicit GetDocSLP_NS(const TStorage& t_storage) : Base(t_storage) {}

  GetDocSLP_NS() = default;

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

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return slp_ ? sdsl::serialize(*slp_, out, child, "slp") : sdsl::serialize_empty_object<TSLP>(out, child, "slp");
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (slp_)
      append(r, "slp", sdsl::size_in_bytes(*slp_));
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    key_slp_ = t_keys[conf::kSLPNS].get<std::string>();
    slp_ = this->template loadItemPtr<TSLP>(key_slp_, t_source, true);
  }

  std::string key_slp_;
  const TSLP* slp_ = nullptr;
};

// Build the bare-SLP cache (kSLPNS) by delegating to dret::construct(grammar::SLP&,
// Config&, datafile). Idempotent: short-circuits if the typed cache file already
// exists. DA is guaranteed by the surrounding RMQ core's construct(), which calls
// EnsureBasicStructures before invoking this through construct(get_doc_policy()).
template <typename TStorage, uint8_t t_width, typename TSLP>
void construct(GetDocSLP_NS<TStorage, t_width, TSLP>& /*unused*/, Config& t_config) {
  using namespace conf;

  const auto key_slp = t_config.keys[kSLPNS].get<std::string>();
  if (sdsl::cache_file_exists<TSLP>(key_slp, t_config))
    return;

  auto event = sdsl::memory_monitor::event(key_slp);
  auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
  TSLP slp;
  dret::construct(slp, t_config, filepath_da);
}

// Document-array lookup over a differential grammar-compressed SLP. The default
// type matches dgcda::DocListIdxDGCDA's default TSLP for typed cache sharing.
template <typename TStorage = GenericStorage, uint8_t t_width = 8, typename TDSLP = DifferentialLightSLP<>>
class GetDocDSLP : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using DSLP = TDSLP;
  using typename Base::size_type;

  explicit GetDocDSLP(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  GetDocDSLP() = default;

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

  uint32_t block_size() const {
    return block_size_;
  }

  float storing_factor() const {
    return storing_factor_;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return dslp_ ? sdsl::serialize(*dslp_, out, child, "dslp")
                 : sdsl::serialize_empty_object<TDSLP>(out, child, "dslp");
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (dslp_)
      collectSizes(r, *dslp_, "dslp_");
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    key_dslp_ =
        std::format("{}-{}_", block_size_, storing_factor_) + t_keys[conf::kDGCDA][conf::kSLP].get<std::string>();
    dslp_ = this->template loadItemPtr<TDSLP>(key_dslp_, t_source, true);
  }

  std::string key_dslp_;
  const TDSLP* dslp_ = nullptr;
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

template <typename TStorage, uint8_t t_width, typename TDSLP>
void construct(GetDocDSLP<TStorage, t_width, TDSLP>& t_get_doc, Config& t_config) {
  using namespace conf;

  const auto key_prefix = std::format("{}-{}_", t_get_doc.block_size(), t_get_doc.storing_factor());
  const auto key_dslp = key_prefix + t_config.keys[kDGCDA][kSLP].get<std::string>();
  if (sdsl::cache_file_exists<TDSLP>(key_dslp, t_config))
    return;

  auto event = sdsl::memory_monitor::event(key_dslp);
  TDSLP dslp;
  dret::construct(dslp, t_config, t_get_doc.block_size(), t_get_doc.storing_factor());
}

}  // namespace dret::rmq
