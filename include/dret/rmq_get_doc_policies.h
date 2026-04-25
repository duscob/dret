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

#pragma once

#include <cstddef>
#include <format>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "config.h"
#include "doc_list_sampled_tree_gcda.h"
#include "index_base.h"
#include "size_report.h"
#include "slp_tools.h"

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

}  // namespace dret::rmq
