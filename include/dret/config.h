//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/9/25.
//

#pragma once

#include "sr-index/alphabet.h"
#include "sr-index/config.h"
#include "sr-index/index_base.h"

namespace dret {

using sri::GenericStorage;

using sri::Config;

using sri::JSON;

using sri::Alphabet;

//~~~~~~~

namespace conf {
constexpr std::string_view kDocEnds = "docEnd";
}

template <uint8_t t_width>
auto createDefaultKeys() {
  auto keys = sri::createDefaultKeys<t_width>();

  keys.update({{conf::kDocEnds, "doc_end"}});

  return keys;
}

}  // namespace dret