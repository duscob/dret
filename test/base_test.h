//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//

#pragma once

#include <cstdint>

#include <gtest/gtest.h>

#include "dret/config.h"


using String = std::string;

template <uint8_t t_width = DRET_DEFAULT_ALPHABET_WIDTH>
class BaseConfigTests : public testing::Test {
 protected:
  template <typename TData>
  void Init(const TData& t_data, sri::SAAlgo t_sa_algo = sri::SAAlgo::SDSL_SE_SAIS) {
    const uint8_t data_width = sizeof(typename TData::value_type) * 8;
    config_ = dret::Config(
        "", std::filesystem::current_path(), t_sa_algo, false, data_width, 1, dret::createDefaultKeys<t_width>());

    auto data_path = sdsl::cache_file_name(key_tmp_input_, config_);
    config_.data_path = data_path;
    std::ofstream out(data_path, std::ios::out | std::ios::binary);
    if (!out) {
      FAIL() << "Cannot open file: " << data_path;
    }
    out.write(reinterpret_cast<const char*>(t_data.data()), sizeof(typename TData::value_type) * t_data.size());
    out.close();
    register_cache_file(key_tmp_input_, config_);
  }

  void TearDown() override {
    sdsl::util::delete_all_files(config_.file_map);
  }

  dret::Config config_;
  std::string key_tmp_input_ = "data";
};
