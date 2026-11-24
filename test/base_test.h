//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//

#pragma once

#include <cstdint>
#include <filesystem>
#include <system_error>

#include <gtest/gtest.h>

#include "dret/config.h"


using String = std::string;

template <uint8_t t_width = DRET_DEFAULT_ALPHABET_WIDTH>
class BaseConfigTests : public testing::Test {
 protected:
  template <typename TData>
  void Init(const TData& t_data, sri::SAAlgo t_sa_algo = sri::SAAlgo::SDSL_SE_SAIS) {
    // Each test gets its own throwaway directory under the OS temp dir.
    // construct()'s internal store_to_cache calls don't go through
    // register_cache_file, so file_map-based cleanup misses them — wiping
    // the whole directory in TearDown collects everything regardless.
    const auto* test_info = testing::UnitTest::GetInstance()->current_test_info();
    std::string subdir = "dret-test-";
    if (test_info) {
      subdir += test_info->test_suite_name();
      subdir += '.';
      subdir += test_info->name();
      subdir += '.';
    }
    subdir += std::to_string(::getpid());
    // Typed-test suite names contain '/', which would create
    // intermediate directories that remove_all of the leaf can't
    // collect. Flatten the path to a single segment.
    for (auto& c : subdir) {
      if (c == '/' || c == '\\') c = '_';
    }
    tmp_dir_ = std::filesystem::temp_directory_path() / subdir;
    std::filesystem::create_directories(tmp_dir_);

    const uint8_t data_width = sizeof(typename TData::value_type) * 8;
    config_ = dret::Config(
        "", tmp_dir_, t_sa_algo, false, data_width, 1, dret::createDefaultKeys<t_width>());

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
    if (!tmp_dir_.empty()) {
      std::error_code ec;
      std::filesystem::remove_all(tmp_dir_, ec);
    }
  }

  dret::Config config_;
  std::filesystem::path tmp_dir_;
  std::string key_tmp_input_ = "data";
};
