//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The suffix array in a differential grammar (DifferentialSLP: Re-Pair over the differences of
// consecutive values, with sampled absolute values every block of positions), as a backend
// (backend/backend.h) that also locates: Locate(i) is SA[i], and the document of position i is the
// number of document ends before SA[i]. Document frequencies over per-document inverse suffix arrays
// read the suffixes from it (doc_freq/oracle.h). The 2020 frequency code's DSA.
//

#pragma once

#include <cstddef>
#include <string>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>

#include "dret/backend/backend.h"
#include "dret/collection.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"
#include "dret/slp/differential_slp.h"
#include "dret/slp/slp_components.h"

namespace dret::backend {

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TDSLP = DifferentialSLP<>,
          typename TBvDocEnds = sdsl::sd_vector<>>
class DiffSA : public FieldsComponent<DiffSA<TStorage, t_width, TDSLP, TBvDocEnds>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<DiffSA, TStorage, t_width>;
  using typename Base::size_type;

  // t_block_size: the spacing of the absolute samples (the cost of a lookup against the space).
  explicit DiffSA(const TStorage& t_storage, uint32_t t_block_size = 256) : Base(t_storage), block_size_(t_block_size) {}

  // As the RMQ cores build a grammar backend: (storage, block size, storing factor); no storing factor.
  DiffSA(const TStorage& t_storage, uint32_t t_block_size, float /*t_storing_factor*/) : DiffSA(t_storage, t_block_size) {}

  DiffSA() = default;

  // SA[t_i].
  std::size_t Locate(std::size_t t_i) const {
    std::size_t value = 0;
    auto report = [&value](auto v) { value = static_cast<std::size_t>(v); };
    ExpandSLP(*dsa_, t_i, t_i + 1, report);
    return value;
  }

  std::size_t operator()(std::size_t t_i) const { return (*rank_)(Locate(t_i)); }

  template <typename TReport>
  void operator()(std::size_t t_b, std::size_t t_e, TReport& t_report) const {
    auto report = [this, &t_report](auto v) { t_report(static_cast<std::size_t>((*rank_)(v))); };
    ExpandSLP(*dsa_, t_b, t_e, report);
  }

  template <typename F>
  void ExpandUntil(std::size_t t_b, std::size_t t_e, F&& t_f) const {
    auto report = [this, &t_f](auto v) { return t_f(static_cast<std::size_t>((*rank_)(v))); };
    ExpandSLPUntil(*dsa_, t_b, t_e, report);
  }

  std::size_t size() const { return dsa_ ? dsa_->SeqSize() : 0; }

  uint32_t block_size() const { return block_size_; }

 protected:
  friend Base;

  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto ends = [](const JSON& k) { return k[conf::kDocEnds].template get<std::string>(); };
    return std::tuple{Item("dsa", t_self.dsa_, [&t_self](const JSON& k) { return DiffSAKey(k, t_self.block_size_); }),
                      Item("doc_end", t_self.doc_ends_, ends),
                      Rank<TBvDocEnds>("doc_end_rank", t_self.rank_, ends)};
  }

 public:
  static std::string DiffSAKey(const JSON& t_keys, uint32_t t_block_size) {
    return PrefixedKey(t_keys, conf::kBlk, KeyName(t_keys, conf::kDiffSA), t_block_size);
  }

 protected:
  const TDSLP* dsa_ = nullptr;
  const TBvDocEnds* doc_ends_ = nullptr;
  const typename TBvDocEnds::rank_1_type* rank_ = nullptr;
  uint32_t block_size_ = 256;
};

template <typename TStorage, uint8_t t_width, typename TDSLP, typename TBvDocEnds>
void construct(DiffSA<TStorage, t_width, TDSLP, TBvDocEnds>& t_backend, Config& t_config) {
  EnsureSA<t_width>(t_config);
  EnsureDocEnds<t_width, TBvDocEnds>(t_config);
  const auto key = DiffSA<TStorage, t_width, TDSLP, TBvDocEnds>::DiffSAKey(t_config.keys, t_backend.block_size());
  if (sdsl::cache_file_exists<TDSLP>(key, t_config)) return;
  auto event = sdsl::memory_monitor::event(key);
  sdsl::int_vector<> sa;
  sdsl::load_from_cache(sa, t_config.keys[conf::kSA].get<std::string>(), t_config);
  TDSLP dsa;
  dsa.Compute(sa, t_backend.block_size());
  sdsl::store_to_cache(dsa, key, t_config, true);
}

}  // namespace dret::backend
