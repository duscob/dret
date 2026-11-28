//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/9/25.
//

#pragma once

#include <cstddef>
#include <vector>

#include "sr-index/sr_idx_generic.h"
#include "sr-index/sr_index.h"

#include "dret/collection.h"
#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/fields.h"
#include "dret/index_base.h"

namespace dret {

template <typename TStorage = GenericStorage, typename TBvDocEnds = sdsl::sd_vector<>, uint8_t t_width = 8>
class GetDocBv;

//~~~~~~~


// Brute force over the r-index or sr-index locate (Brute-RI / Brute-SRI): the text position of every
// occurrence, mapped to its document by the document ends (GetDocBv). A pattern core: the locate index
// finds the occurrences itself, so the index holds no r-index count (doc_list_index.h). Serialized as
// the locate index, then the document ends; no size report (the benchmark falls back to the serialized
// size).
template <typename TStorage, typename TAlphabet, typename TLocateIdx, typename TGetDoc>
class LocateCore
    : public FieldsComponent<LocateCore<TStorage, TAlphabet, TLocateIdx, TGetDoc>, TStorage, TAlphabet::int_width> {
 public:
  using Base = FieldsComponent<LocateCore, TStorage, TAlphabet::int_width>;
  using typename Base::size_type;
  using TPattern = typename TAlphabet::string_type;

  // Every occurrence goes to the sink once.
  static constexpr bool kFrequencies = true;

  explicit LocateCore(const TStorage& t_storage) : Base(t_storage), locate_idx_(t_storage), get_doc_(t_storage) {}

  LocateCore(const TStorage& t_storage, const TLocateIdx& t_locate_idx)
      : Base(t_storage), locate_idx_(t_locate_idx), get_doc_(t_storage) {}

  LocateCore(const TStorage& t_storage, const TLocateIdx& t_locate_idx, const TGetDoc& t_get_doc)
      : Base(t_storage), locate_idx_(t_locate_idx), get_doc_(t_get_doc) {}

  LocateCore() = default;

  // The document of every occurrence Locate yields.
  template <typename TSink>
  void Query(const TPattern& t_pattern, TSink& t_sink) const {
    for (const auto& item : locate_idx_.Locate(t_pattern)) t_sink.Position(get_doc_(item));
  }

  std::size_t n_doc() const { return get_doc_.n_doc(); }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h). Neither is sized:
  // the benchmark then reports the serialized size (Brute-RI / Brute-SRI as published).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    return std::tuple{SriPart<SizeAs::kNone>("locate_idx", t_self.locate_idx_),
                      Part<SizeAs::kNone>("get_doc", t_self.get_doc_)};
  }

  TLocateIdx locate_idx_;
  TGetDoc get_doc_;
};

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TLocateIdx = sri::SrIdxGeneric<sri::SrIndexValidArea<TStorage, TAlphabet>, 16>,
          typename TGetDoc = GetDocBv<TStorage, sdsl::sd_vector<>, TAlphabet::int_width>>
using DocListIdxBrute =
    DocListIdx<TStorage, TAlphabet, NoCountIdx, LocateCore<TStorage, TAlphabet, TLocateIdx, TGetDoc>>;

//~~~~~~~


template <typename TStorage, typename TAlphabet, typename TLocateIdx, typename TGetDoc>
void construct(LocateCore<TStorage, TAlphabet, TLocateIdx, TGetDoc>& t_core, Config& t_config) {
  EnsureText<TAlphabet::int_width>(t_config);

  TLocateIdx locate_index(t_core.storage());
  construct(locate_index, t_config.data_path, t_config);

  TGetDoc get_doc(t_core.storage());
  construct(get_doc, t_config);
}

// Overload for locate indexes that take the subsample rate as a runtime
// constructor argument (the sri::SrIndex family) rather than as a template
// parameter (sri::SrIdxGeneric, which is only a thin wrapper forwarding a
// compile-time rate to this same constructor).
//
// Construct mode needs this because the spec's sampling_size is a list of
// runtime values; binding the rate at compile time would restrict it to a
// fixed set. The rate selects the "<rate>_" cache-key prefix
// (sri::SrIndex::setupKeyNames), so it must match the rate the query side
// later loads with, or the load fails with a missing "<rate>_bwt_run_*" key.
//
// The index must already carry the same rate: construct(index, config, rate)
// ends by load()ing the index's own locate index, and a default-constructed one
// has an empty key prefix.
template <typename TStorage, typename TAlphabet, typename TLocateIdx, typename TGetDoc>
void construct(LocateCore<TStorage, TAlphabet, TLocateIdx, TGetDoc>& t_core,
               Config& t_config,
               std::size_t t_subsample_rate) {
  EnsureText<TAlphabet::int_width>(t_config);

  TLocateIdx locate_index(t_core.storage(), t_subsample_rate);
  construct(locate_index, t_config.data_path, t_config);

  TGetDoc get_doc(t_core.storage());
  construct(get_doc, t_config);
}

//~~~~~~~


template <typename TStorage, typename TBvDocEnds, uint8_t t_width>
class GetDocBv : public FieldsComponent<GetDocBv<TStorage, TBvDocEnds, t_width>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<GetDocBv, TStorage, t_width>;
  using typename Base::size_type;

  GetDocBv(const TStorage& t_storage) : Base(t_storage) {}

  GetDocBv() = default;

  uint32_t operator()(uint32_t t_position) const {
    return get_doc_(t_position);
  }

  // The number of documents: the ones of the document-ends bitvector.
  std::size_t n_doc() const { return n_doc_; }

 protected:
  friend Base;

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    auto key = [](const JSON& k) { return k[conf::kDocEnds].template get<std::string>(); };
    return std::tuple{Item("doc_end", t_self.doc_ends_, key), Rank<TBvDocEnds>("doc_end_rank", t_self.rank_, key)};
  }

  void OnLoad() {
    n_doc_ = (*rank_)(doc_ends_->size());
    get_doc_ = [rank = rank_](uint32_t t_position) { return (*rank)(t_position); };
  }

  //~~~~~~~

  const TBvDocEnds* doc_ends_ = nullptr;
  const typename TBvDocEnds::rank_1_type* rank_ = nullptr;
  std::function<uint32_t(uint32_t)> get_doc_;
  std::size_t n_doc_ = 0;
};

//~~~~~~~


template <typename TStorage, typename TBvDocEnds, uint8_t t_width>
void construct(GetDocBv<TStorage, TBvDocEnds, t_width>& t_index, Config& t_config) {
  EnsureDocEnds<t_width, TBvDocEnds>(t_config);

  t_index.load(t_config);
}


}  // namespace dret