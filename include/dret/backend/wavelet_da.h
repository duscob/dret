//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The document array in a wavelet tree, as a backend (backend/backend.h) that also counts: Rank(i, d),
// the occurrences of d in DA[0, i). Document frequencies use it (doc_freq/oracle.h): a document's
// frequency in [l, r] is Rank(r + 1, d) - Rank(l, d). The 2020 frequency code's WT-DA, sdsl::wt_ap.
//

#pragma once

#include <cstddef>
#include <string>

#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "dret/backend/backend.h"
#include "dret/collection.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret::backend {

template <typename TStorage = GenericStorage, uint8_t t_width = 8, typename TWt = sdsl::wt_ap<>>
class WaveletDA : public FieldsComponent<WaveletDA<TStorage, t_width, TWt>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<WaveletDA, TStorage, t_width>;
  using typename Base::size_type;

  explicit WaveletDA(const TStorage& t_storage) : Base(t_storage) {}

  WaveletDA() = default;

  std::size_t operator()(std::size_t t_i) const { return static_cast<std::size_t>((*wt_)[t_i]); }

  template <typename TReport>
  void operator()(std::size_t t_b, std::size_t t_e, TReport& t_report) const {
    for (std::size_t i = t_b; i < t_e; ++i) t_report(static_cast<std::size_t>((*wt_)[i]));
  }

  template <typename F>
  void ExpandUntil(std::size_t t_b, std::size_t t_e, F&& t_f) const {
    for (std::size_t i = t_b; i < t_e; ++i)
      if (!t_f(static_cast<std::size_t>((*wt_)[i]))) return;
  }

  // The occurrences of t_doc in DA[0, t_i).
  std::size_t Rank(std::size_t t_i, std::size_t t_doc) const { return wt_->rank(t_i, t_doc); }

  std::size_t size() const { return wt_ ? wt_->size() : 0; }

 protected:
  friend Base;

  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    return std::tuple{Item("wt", t_self.wt_, [](const JSON& k) { return k[conf::kDaWt].template get<std::string>(); })};
  }

  const TWt* wt_ = nullptr;
};

template <typename TStorage, uint8_t t_width, typename TWt>
void construct(WaveletDA<TStorage, t_width, TWt>& /*unused*/, Config& t_config) {
  EnsureDA<t_width>(t_config);
  const auto key = t_config.keys[conf::kDaWt].get<std::string>();
  if (sdsl::cache_file_exists<TWt>(key, t_config)) return;
  auto event = sdsl::memory_monitor::event(key);
  // The DA's file is a serialized int_vector (num_bytes 0).
  TWt wt;
  sdsl::construct(wt, sdsl::cache_file_name<sdsl::int_vector<>>(t_config.keys[conf::kDA].get<std::string>(), t_config), 0);
  sdsl::store_to_cache(wt, key, t_config, true);
}

}  // namespace dret::backend
