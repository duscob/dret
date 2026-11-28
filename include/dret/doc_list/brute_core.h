//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// Brute force over a backend: the document of every position of [sp, ep), read by the backend
// (backend/backend.h), into the sink, which keeps each document once. No precomputation. Over the
// plain DA, SA-Phi or the RLCSA it is Brute-DA / Brute-SAPhi / Brute-RLCSA (doc_list_brute_backend.h);
// over a GCDA grammar without its lists, Brute-GCDA (doc_list_slp.h).
//
// A backend that finds the range itself (the RLCSA) does so for the index (Count), which then holds no
// r-index (doc_list_index.h). The RLCSA must count for itself anyway, since its suffix array orders
// suffixes differently around the separators (see factories/rmq.h).
//
// Serialized as the backend, then the number of documents.
//

#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <variant>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "dret/backend/backend.h"
#include "dret/collection.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret {

template <typename TStorage, uint8_t t_width, typename TBackend>
class BruteCore : public FieldsComponent<BruteCore<TStorage, t_width, TBackend>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<BruteCore, TStorage, t_width>;
  using typename Base::size_type;
  using Backend = TBackend;

  // The backend's parameters follow the storage (e.g. a grammar's block size and storing factor).
  template <typename... TArgs>
    requires std::constructible_from<TBackend, const TStorage&, TArgs...>
  explicit BruteCore(const TStorage& t_storage, TArgs&&... t_backend_args)
      : Base(t_storage), backend_(t_storage, std::forward<TArgs>(t_backend_args)...) {}

  BruteCore() = default;

  // A reference to the sink, not the sink: some expansions take their report by value.
  template <typename TSink>
  void Query(std::size_t t_sp, std::size_t t_ep, std::size_t /*t_m*/, TSink& t_sink) const {
    auto position = [&t_sink](auto d) { t_sink.Position(static_cast<std::size_t>(d)); };
    backend_(t_sp, t_ep, position);
  }

  template <typename TPattern>
    requires ::dret::backend::CountsPattern<TBackend, TPattern>
  std::pair<std::size_t, std::size_t> Count(const TPattern& t_pattern) const {
    return backend_.Count(t_pattern);
  }

  std::size_t n_doc() const { return n_doc_; }

  const TBackend& backend() const { return backend_; }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    return std::tuple{Part("get_doc", t_self.backend_), NDoc(t_self.n_doc_)};
  }

  TBackend backend_;
  std::size_t n_doc_ = 0;
};

// The collection's structures (EnsureCollection: text, SA, document ends, DA and the number of
// documents), then the backend's own. A copy of the backend, as TBackend need not be
// default-constructible (not when TStorage is std::reference_wrapper).
template <typename TStorage, uint8_t t_width, typename TBackend>
void construct(BruteCore<TStorage, t_width, TBackend>& t_core, Config& t_config) {
  EnsureCollection<t_width>(t_config);

  auto backend = t_core.backend();
  construct(backend, t_config);
}

}  // namespace dret
