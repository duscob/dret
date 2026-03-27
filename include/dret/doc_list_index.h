//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/2/21.
//
#pragma once

#include <functional>
#include <memory>
#include <string>

#include "config.h"
#include "index_base.h"

namespace dret {

template <typename TSequence = Alphabet<>::string_type>
class DocListIndex {
 public:
  using TPattern = TSequence;
  using TDocId = std::size_t;

  virtual ~DocListIndex() = default;

  virtual void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const = 0;
};

//~~~~~~~


template <typename TStorage = GenericStorage, typename TAlphabet = Alphabet<>>
class DocListIndexExtStorage : public DocListIndex<typename TAlphabet::string_type>,
                               public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  DocListIndexExtStorage(const TStorage& t_storage) : IndexBaseWithExternalStorage<TStorage>(t_storage) {}

  DocListIndexExtStorage() = default;
};

//~~~~~~~


template <typename TLocate, typename TGetDoc, typename TSequence = Alphabet<>::string_type>
class DocListIndexBrute : public DocListIndex<TSequence> {
 public:
  using Base = DocListIndex<TSequence>;
  using typename Base::TDocId;
  using typename Base::TPattern;

  DocListIndexBrute(const TLocate& t_locate, const TGetDoc& t_get_doc) : locate_{t_locate}, get_doc_{t_get_doc} {}

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto occurrences = locate_(t_pattern);

    for (const auto& item : occurrences) {
      t_report(get_doc_(item));
    }
  }

 private:
  TLocate locate_;
  TGetDoc get_doc_;
};

//~~~~~~~


template <typename TComputeSARange, typename TComputeDocs, typename TSequence = Alphabet<>::string_type>
class DocListIndexBasicScheme : public DocListIndex<TSequence> {
 public:
  using Base = DocListIndex<TSequence>;
  using typename Base::TDocId;
  using typename Base::TPattern;

  DocListIndexBasicScheme(const TComputeSARange& t_csa, const TComputeDocs& t_compute_docs)
      : compute_sa_range_{t_csa}, compute_docs_{t_compute_docs} {}

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = compute_sa_range_(t_pattern);

    compute_docs_(sp, ep, t_report);
  }

 private:
  TComputeSARange compute_sa_range_;
  TComputeDocs compute_docs_;
};

}  // namespace dret
