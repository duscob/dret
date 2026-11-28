//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The cells of a sweep: one per index the spec asks for, the same in both modes. A family's source
// file (factories/<family>.cpp) turns its sweep block into cells: the name (FAMILY[k=v|...]), how to
// make the index, and what construct mode needs besides (the cache files --rebuild deletes, the
// counters). Query mode builds what is missing and loads; construct mode times construct(). So a
// family states its cells, their names and their types once, and both modes read them.
//
// Each family is compiled in a source file of its own (RMQ one per core, PDL one per codec): in one
// translation unit holding every family, GCC's inlining decisions in one family's query code
// depended on how much code the others added (the unit-growth budget), so a change anywhere could
// move a family's time by a few percent (seen 2026-11 on Brute-DA, CountRI and PDL, with sdsl's
// rank_support_v::rank inlined in one build and not in the next). Split, a family's object code
// changes only when its own code, or a header it includes, does.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/io.hpp>

#include "sr-index/index_base.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_freq/frequencies.h"
#include "dret/doc_list/doc_index.h"

#include "../spec.h"

namespace bench::factories {

// Every index lives in a storage it shares with the other indexes of the run (query mode) or that
// its benchmark owns (construct mode).
using Storage = std::reference_wrapper<sri::GenericStorage>;

// t_report(d, f) for each document d of a pattern, with its frequency f.
using FreqReport = std::function<void(std::size_t, std::size_t)>;
using Frequencies = std::function<void(const std::string&, const FreqReport&)>;

// An index, made but neither built nor loaded, behind its type-erased interface, with the calls
// that need its type.
struct CellIndex {
  std::shared_ptr<dret::DocListIndex<>> index;
  std::function<void(dret::Config&)> construct;  // builds the cache files the index needs, if missing
  std::function<void(dret::Config&)> load;
  std::function<std::size_t()> serialized_size;
  Frequencies frequencies;  // empty when the index does not compute document frequencies
};

// The index computes document frequencies (dret/doc_freq/frequencies.h).
template <typename TIndex>
concept ComputesFrequencies = requires {
  typename TIndex::Core;
  requires dret::freq::FrequencyIndex<TIndex>;
};

using MakeIndex = std::function<CellIndex(Storage)>;

// The frequency query of t_idx, or nothing when it does not compute frequencies. Every cell index
// sets its frequencies from this.
template <typename TIndex>
Frequencies FrequenciesOf(const std::shared_ptr<TIndex>& t_idx) {
  if constexpr (ComputesFrequencies<TIndex>) {
    return [t_idx](const std::string& t_pattern, const FreqReport& t_report) {
      dret::freq::Frequencies(*t_idx, t_pattern, t_report);
    };
  } else {
    return {};
  }
}

// The cell index of TIndex as t_new(storage) creates it.
template <typename TIndex, typename TNew>
MakeIndex MakeBy(TNew t_new) {
  return [t_new](Storage t_storage) {
    std::shared_ptr<TIndex> idx = t_new(t_storage);
    return CellIndex{idx,
                     [idx](dret::Config& t_config) { construct(*idx, t_config); },
                     [idx](dret::Config& t_config) { idx->load(t_config); },
                     [idx] { return static_cast<std::size_t>(sdsl::size_in_bytes(*idx)); },
                     FrequenciesOf(idx)};
  };
}

// The cell index of TIndex(storage, t_args...).
template <typename TIndex, typename... TArgs>
MakeIndex MakeOf(TArgs... t_args) {
  return MakeBy<TIndex>([=](Storage t_storage) { return std::make_shared<TIndex>(t_storage, t_args...); });
}

struct Cell {
  std::string name;
  MakeIndex make;
  // Query mode builds a missing index before loading it, except the brute baselines over the r-index
  // and the sr-index: their construct() rebuilds the locate index from the text.
  bool query_constructs = true;
  // Construct mode with workload.rebuild: the cache-key prefixes of the files this cell owns, deleted
  // before each timed construct(). Empty: no rebuild (the files are shared or managed elsewhere).
  std::vector<std::string> rebuild_prefixes;
  // Construct-mode counters (bs, sf, s), when not 0.
  std::uint32_t block_size = 0;
  float storing_factor = 0;
  std::size_t sampling = 0;
};

// One function per sweep block, defined in the family's source file.
std::vector<Cell> Cells(const spec::BruteSweep& t_sweep);        // brute.cpp
std::vector<Cell> Cells(const spec::GCDASweep& t_sweep);         // gcda.cpp
std::vector<Cell> Cells(const spec::GCDABackendSweep& t_sweep);  // gcda_backend.cpp
std::vector<Cell> Cells(const spec::DGCDASweep& t_sweep);        // dgcda.cpp
std::vector<Cell> Cells(const spec::SLPNSSweep& t_sweep);        // slp_ns.cpp
std::vector<Cell> Cells(const spec::RMQSweep& t_sweep);          // rmq.cpp (+ rmq_<core>.cpp)
std::vector<Cell> Cells(const spec::PDLSweep& t_sweep);          // pdl.cpp (+ pdl_<codec>.cpp)
std::vector<Cell> Cells(const spec::RMQFreqSweep& t_sweep);      // rmq_freq.cpp
std::vector<Cell> Cells(const spec::GCDAOccsSweep& t_sweep);     // gcda_occs.cpp

// Occurrences of a pattern (the width of its range in the r-index), for the per-pattern counters.
// Loads the r-index on its first call, so it needs the same cache files as Brute-RI. (brute.cpp)
std::function<std::size_t(const std::string&)> OccurrenceCounter(Storage t_storage, const dret::Config& t_config);

}  // namespace bench::factories
