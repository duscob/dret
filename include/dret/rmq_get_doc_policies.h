//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// GetDoc policies for RMQ document-listing indexes. Each policy implements
// single-position and range-based document-array lookup via a common interface:
//
//   std::size_t operator()(std::size_t i) const     — doc at SA position i
//   void operator()(std::size_t b, std::size_t e, TReport&) const  — half-open [b,e)
//   std::size_t size() const                        — total SA length
//
// GetDocDA: plain sdsl::int_vector<> document array (no grammar compression).

#pragma once

#include <cstddef>
#include <string>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "config.h"
#include "index_base.h"
#include "size_report.h"

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

}  // namespace dret::rmq
