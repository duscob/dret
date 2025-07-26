//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <any>
#include <map>
#include <string>

#include "sr-index/config.h"

namespace dret {

using GenericStorage = std::map<std::string, std::any>;

using Config = sri::Config;

}  // namespace dret
