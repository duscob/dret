//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The brute-force cells: over the r-index (Brute-RI), the sr-index (Brute-SRI, one per subsample
// rate), the plain DA (Brute-DA), the RLCSA (Brute-RLCSA), and the r-index count alone (CountRI).
//

#include <memory>
#include <string>
#include <vector>

#include "brute.h"
#include "cells.h"
#include "names.h"

namespace bench::factories {

namespace {

// Brute-SRI: the subsample rate is a runtime value that selects the "<rate>_" cache-key prefix, so the
// index carries it in its locate index and its construct() takes it too.
MakeIndex MakeBruteSr(std::size_t t_rate) {
  using TIndex = brute::IdxSr<Storage>;
  return [t_rate](Storage t_storage) {
    auto idx = std::make_shared<TIndex>(t_storage, sri::SrIndexValidArea<Storage>(t_storage, t_rate));
    return CellIndex{idx,
                     [idx, t_rate](dret::Config& t_config) { construct(*idx, t_config, t_rate); },
                     [idx](dret::Config& t_config) { idx->load(t_config); },
                     [idx] { return static_cast<std::size_t>(sdsl::size_in_bytes(*idx)); }};
  };
}

}  // namespace

std::vector<Cell> Cells(const spec::BruteSweep& t_sweep) {
  using Kind = spec::BruteSweep::Kind;
  std::vector<Cell> cells;
  for (auto kind : t_sweep.kind) {
    switch (kind) {
      case Kind::RIndex:
        cells.push_back({.name = "BruteRI", .make = MakeOf<brute::IdxR<Storage>>(), .query_constructs = false});
        break;
      case Kind::SrIndex:
        for (auto s : t_sweep.sampling_size) {
          Cell cell{.name = names::KeyedName("BruteSRI", {{"sampling", names::IntStr(s)}}),
                    .make = MakeBruteSr(s),
                    .query_constructs = false};
          cell.sampling = s;
          cells.push_back(std::move(cell));
        }
        break;
      case Kind::DA:
        cells.push_back({.name = "BruteDA", .make = MakeOf<brute::IdxDA<Storage>>()});
        break;
      case Kind::RLCSA:
        cells.push_back({.name = "BruteRLCSA", .make = MakeOf<brute::IdxRLCSA<Storage>>()});
        break;
      case Kind::Count:
        cells.push_back({.name = "CountRI", .make = MakeOf<brute::IdxCount<Storage>>()});
        break;
    }
  }
  return cells;
}

std::function<std::size_t(const std::string&)> OccurrenceCounter(Storage t_storage, const dret::Config& t_config) {
  std::shared_ptr<sri::RIndex<Storage>> r_index;
  return [t_storage, t_config, r_index](const std::string& t_pattern) mutable -> std::size_t {
    if (!r_index) {
      r_index = std::make_shared<sri::RIndex<Storage>>(t_storage);
      r_index->load(t_config);
    }
    const auto [start, end] = r_index->Count(t_pattern);
    return end > start ? end - start : 0;
  };
}

}  // namespace bench::factories
