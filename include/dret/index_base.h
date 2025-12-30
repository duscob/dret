//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <any>
#include <functional>
#include <map>
#include <string>
#include <variant>

#include <sdsl/io.hpp>

#include "config.h"

namespace dret {


template <typename TIndex>
void construct(TIndex& t_index, const std::string& t_data_path, Config& t_config) {
  constructItems(t_index, t_config);

  // t_index.load(t_config);
}

}  // namespace dret
