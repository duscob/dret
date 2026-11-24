//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/2/21.
//
#pragma once

#include <functional>
#include <memory>
#include <string>

#include "dret/config.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret {

template <typename TSequence = Alphabet<>::string_type>
class DocListIndex {
 public:
  using TPattern = TSequence;
  using TDocId = std::size_t;

  virtual ~DocListIndex() = default;

  virtual void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const = 0;

  // Per-field size breakdown for benchmark reporting. Default empty; concrete
  // indices override to report their internal structure sizes.
  virtual SizeReport GetSizeReport() const { return {}; }
};

//~~~~~~~


template <typename TStorage = GenericStorage, typename TAlphabet = Alphabet<>>
class DocListIndexExtStorage : public DocListIndex<typename TAlphabet::string_type>,
                               public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  DocListIndexExtStorage(const TStorage& t_storage) : IndexBaseWithExternalStorage<TStorage>(t_storage) {}

  DocListIndexExtStorage() = default;
};

}  // namespace dret
