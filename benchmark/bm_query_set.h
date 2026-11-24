//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/8/26.
//

#pragma once

#include <string>

#include "bm_query.h"

template <typename TFactoryConfig>
struct IndexConfig {
  std::string name;
  TFactoryConfig factory_config;
  bool has_sampling = false;
};

//~~~~~~~


auto RegisterAllQueryBenchmarks = [](auto t_factory,
                                     auto t_n,
                                     const auto& t_idx_configs,
                                     const auto& t_patterns,
                                     auto t_bm_config,
                                     auto t_min_s,
                                     auto t_max_s,
                                     auto t_multiplier_s = 2) {
  for (auto idx_config : t_idx_configs) {
    // Index builder
    auto make_index = [t_factory, idx_config](const auto& tt_state) {
      auto fac_config = idx_config.factory_config;
      if (idx_config.has_sampling) {
        fac_config.sampling_size = tt_state.range(0);
      }

      auto idx = t_factory(fac_config);

      return idx;
    };

    // Benchmarks
    auto bms =
        RegisterQueryBenchmarks(idx_config.name, make_index, t_patterns, t_n, t_bm_config, idx_config.has_sampling);
    // Index with sampling
    if (idx_config.has_sampling) {
      for (auto& bm : bms) {
        bm.second->RangeMultiplier(t_multiplier_s)->Range(t_min_s, t_max_s);
      }
    }
  }
};
