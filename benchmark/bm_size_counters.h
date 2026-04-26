//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/24/26.
//

#pragma once

#include <string>
#include <string_view>

#include <benchmark/benchmark.h>

#include "dret/size_report.h"

inline void appendCounters(benchmark::State& state,
                           const dret::SizeReport& report,
                           std::string_view suffix = "_bytes") {
  for (const auto& [name, bytes] : report) {
    state.counters[name + std::string(suffix)] = static_cast<double>(bytes);
  }
}
