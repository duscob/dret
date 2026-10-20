//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/19/26.
//
// Brute force over a backend: find the suffix-array range [sp, ep) of the pattern, read the document
// of every position in it from the backend, and report each document once (distinct_docs.h). No
// precomputation; the backend is any raw-range get-doc policy of the RMQ cores and PDL (the plain DA,
// SA-Phi, a GCDA grammar, the RLCSA), expanded with its range operator.
//
// A backend that finds the range itself (the RLCSA) replaces the r-index, as in PDL: count_idx is then
// neither loaded, serialized nor counted in the size. The RLCSA must count for itself anyway, since its
// suffix array orders suffixes differently around the separators (see factories/rmq.h).
//

#pragma once

#include <cstddef>
#include <functional>
#include <string>

#include <sdsl/int_vector.hpp>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/doc_list/distinct_docs.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/index_base.h"
#include "dret/rmq/rmq_get_doc_policies.h"

namespace dret {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TGetDoc = rmq::GetDocDA<TStorage, TAlphabet::int_width>>
class DocListIdxBruteBackend : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using StorageBase = DocListIndexExtStorage<TStorage, TAlphabet>;
  using StorageBaseImpl = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using TPattern = typename TAlphabet::string_type;
  using TDocId = std::size_t;
  using size_type = std::size_t;

  static constexpr bool kBackendCounts = requires(const TGetDoc& g, const TPattern& p) { g.Count(p); };

  explicit DocListIdxBruteBackend(const TStorage& t_storage)
      : StorageBase(t_storage), count_idx_(t_storage), get_doc_(t_storage) {}

  DocListIdxBruteBackend() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = Count(t_pattern);
    if (sp >= ep) return;

    DistinctDocs seen(n_doc_);
    auto add = [&seen, &t_report](std::size_t d) {
      if (seen.Insert(d)) t_report(static_cast<TDocId>(d));
    };
    get_doc_(sp, ep, add);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    if constexpr (!kBackendCounts) written += count_idx_.serialize(out, child, "count_idx");
    written += get_doc_.serialize(out, child, "get_doc");
    sdsl::int_vector<64> n_doc_vec(1, n_doc_);
    written += sdsl::serialize(n_doc_vec, out, child, "n_doc");
    return written;
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    if constexpr (!kBackendCounts) append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    for (const auto& [k, v] : get_doc_.GetSizeReport()) append(r, k, v);
    return r;
  }

  const TCountIdx& count_idx() const { return count_idx_; }

  const TGetDoc& get_doc() const { return get_doc_; }

 protected:
  std::pair<std::size_t, std::size_t> Count(const TPattern& t_pattern) const {
    if constexpr (kBackendCounts) return get_doc_.Count(t_pattern);
    else return count_idx_.Count(t_pattern);
  }

  void loadInner(typename StorageBaseImpl::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit(
        [this](auto&& tt_source) {
          if constexpr (!kBackendCounts) count_idx_.load(tt_source.get());
          get_doc_.load(tt_source.get());
        },
        t_source);

    const auto* n_doc_vec =
        this->template loadItemPtr<sdsl::int_vector<64>>(t_keys[kRmqNDoc].get<std::string>(), t_source, true);
    n_doc_ = static_cast<std::size_t>((*n_doc_vec)[0]);
  }

  TCountIdx count_idx_;
  TGetDoc get_doc_;
  std::size_t n_doc_ = 0;
};

//~~~~~~~


// Text, SA, document ends, DA and the number of documents as the RMQ cores build them, the backend's
// own structures and, unless the backend counts, the r-index.
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TGetDoc>
void construct(DocListIdxBruteBackend<TStorage, TAlphabet, TCountIdx, TGetDoc>& t_index, Config& t_config) {
  rmq::internal::EnsureBasicStructures<TAlphabet::int_width, sdsl::sd_vector<>>(t_config);

  auto get_doc = t_index.get_doc();
  construct(get_doc, t_config);

  if constexpr (!DocListIdxBruteBackend<TStorage, TAlphabet, TCountIdx, TGetDoc>::kBackendCounts) {
    auto count_idx = t_index.count_idx();
    construct(count_idx, t_config.data_path, t_config);
  }

  t_index.load(t_config);
}

}  // namespace dret
