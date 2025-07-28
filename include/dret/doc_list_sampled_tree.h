//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <filesystem>

#include "sdsl/io.hpp"

#include "construct_base.h"
#include "doc_list_index.h"
#include "index_base.h"

namespace dret {

template <typename TStorage = GenericStorage>
class DLSampledTreeScheme : public DocListIndex {
 public:
  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {}
};

//~~~~~~~


template <typename TStorage = GenericStorage>
class GCDA : public DLSampledTreeScheme<TStorage> {
 public:
  using Base = DLSampledTreeScheme<TStorage>;
  using typename Base::TDocId;
  using typename Base::TPattern;

  GCDA() = default;

  const uint32_t& block_size() const { return block_size_; }

  const float& storing_factor() const { return storing_factor_; }

 protected:
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

//~~~~~~~


template <typename TStorage = GenericStorage>
void constructItems(DLSampledTreeScheme<TStorage>& t_index, Config& t_config) {
  if (!sdsl::cache_file_exists(sdsl::conf::KEY_SA, t_config)) {
    auto event = sdsl::memory_monitor::event("SA");
    sdsl::construct_sa<8>(t_config);
  }

  if (!cache_file_exists(dret::conf::KEY_DOC_END, t_config)) {
    auto event = sdsl::memory_monitor::event("DocEnds");
    ConstructDocEnd(t_config);
  }

  if (!cache_file_exists(dret::conf::KEY_DA, t_config)) {
    auto event = sdsl::memory_monitor::event("DA");
    ConstructDocArray(t_config);
  }

  if (std::string file_da = cache_file_name(conf::KEY_DA_RAW, t_config);
      !std::filesystem::exists(file_da + ".R") && REPAIR_EXE) {
    auto event = sdsl::memory_monitor::event("Repair DA");
    std::string cmd = REPAIR_EXE + (" " + file_da);
    std::system(cmd.c_str());
  }
}

//~~~~~~~

}  // namespace dret
