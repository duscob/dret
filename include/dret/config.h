//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/9/25.
//

#pragma once

#include "sr-index/alphabet.h"
#include "sr-index/config.h"
#include "sr-index/index_base.h"

namespace dret {

using sri::GenericStorage;

using sri::JSON;

using sri::Alphabet;

//~~~~~~~

namespace conf {
using namespace sri::conf;
constexpr std::string_view kDocEnds = "docEnd";
constexpr std::string_view kDA = "da";
constexpr std::string_view kGCDA = "gcda";
constexpr std::string_view kSLP = "slp";
constexpr std::string_view kDocs = "docs";
}  // namespace conf

template <uint8_t t_width = DRET_DEFAULT_ALPHABET_WIDTH>
struct Keys {
  Keys() {
    keys = sri::createDefaultKeys<t_width>();

    keys.update({
        {conf::kDocEnds, "doc_end"},
        {conf::kDA, "da"},
        {
            conf::kGCDA,
            {
                {conf::kSLP, "gcda_slp"},
                {conf::kDocs, "gcda_docs"},
            },
        },
    });
  }

  nlohmann::json keys;
};

const Keys<> kDefaultKeys;

template <uint8_t t_width>
auto createDefaultKeys() {
  return Keys<t_width>().keys;
}

struct Config : public sri::Config {
  uint64_t data_delim = 0;

  Config() = default;

  Config(const std::filesystem::path& t_data_path,
         const std::filesystem::path& t_output_dir,
         sri::SAAlgo t_sa_algo,
         bool t_delete_files = false,
         uint8_t t_data_width = 8,
         uint64_t t_data_delim = 0,
         JSON t_keys = kDefaultKeys.keys)
      : sri::Config(t_data_path, t_output_dir, t_sa_algo, t_delete_files, t_data_width, std::move(t_keys)),
        data_delim(t_data_delim) {}
};

}  // namespace dret
