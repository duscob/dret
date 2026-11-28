//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// The plain document array as a backend (see backend/backend.h for the contract).
//

#pragma once

#include <cstddef>
#include <string>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "dret/backend/backend.h"
#include "dret/collection.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret::backend {

// Document-array lookup over a plain sdsl::int_vector<>.
// Loads the cache entry at conf::kDA, which EnsureDA builds (dret/collection.h).
template <typename TStorage = GenericStorage, uint8_t t_width = 8>
class PlainDA : public FieldsComponent<PlainDA<TStorage, t_width>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<PlainDA, TStorage, t_width>;
  using typename Base::size_type;

  explicit PlainDA(const TStorage& t_storage) : Base(t_storage) {}

  PlainDA() = default;

  std::size_t operator()(std::size_t i) const {
    return static_cast<std::size_t>((*da_)[i]);
  }

  // Half-open range expansion: calls r(doc) for each position in [b, e).
  template <typename TReport>
  void operator()(std::size_t b, std::size_t e, TReport& r) const {
    for (std::size_t i = b; i < e; ++i)
      r(static_cast<std::size_t>((*da_)[i]));
  }

  // In-order expansion of [b, e) that stops once f returns false (see
  // ExpandSLPUntil in slp_tools.h for why the RMQ cores need it).
  template <typename F>
  void ExpandUntil(std::size_t b, std::size_t e, F&& f) const {
    for (std::size_t i = b; i < e; ++i)
      if (!f(static_cast<std::size_t>((*da_)[i])))
        return;
  }

  // Total number of SA positions. Used by IlcpLikeCore for the last-run sentinel.
  std::size_t size() const {
    return da_ ? da_->size() : 0;
  }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    return std::tuple{Item("da", t_self.da_, [](const JSON& k) { return k[conf::kDA].template get<std::string>(); })};
  }

  const sdsl::int_vector<>* da_ = nullptr;
};

// Builds the DA if missing (dret/collection.h); the indexes over it have built it already.
template <typename TStorage, uint8_t t_width>
void construct(PlainDA<TStorage, t_width>& /*unused*/, Config& t_config) {
  EnsureDA<t_width>(t_config);
}

}  // namespace dret::backend
