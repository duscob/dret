//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode entry points of the factories, each defined in a source file of its own
// (factories/<family>.cpp; RMQ one per core, PDL one per codec). The query factory (factory.h)
// sees only these declarations, so the code an index runs in its timed Search is compiled in that
// family's source file and nowhere else.
//
// Why: in one translation unit holding every family, GCC's inlining decisions in one family's
// query code depended on how much code the other families added (the unit-growth budget), so a
// change anywhere could move a family's time by a few percent (seen 2026-11 on Brute-DA, CountRI
// and PDL, with sdsl's rank_support_v::rank inlined in one build and not in the next). Split, a
// family's object code changes only when its own code, or a header it includes, does.
//
// Each function builds the index's missing cache files (construct) and loads it, as the template
// Make() functions of the family headers do; the brute baselines over the r-index and the sr-index
// only load.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "sr-index/index_base.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"

#include "../axes.h"

namespace bench::factories::query {

using Storage = std::reference_wrapper<sri::GenericStorage>;

// The index behind its type-erased interface, and its serialized size.
using Built = std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>;

// brute.cpp
Built BruteR(Storage t_storage, dret::Config& t_config);
Built BruteSr(Storage t_storage, dret::Config& t_config, std::size_t t_sampling);
Built BruteDA(Storage t_storage, dret::Config& t_config);
Built BruteRLCSA(Storage t_storage, dret::Config& t_config);
Built CountRI(Storage t_storage, dret::Config& t_config);
// Occurrences of a pattern (the width of its range in the r-index), for the per-pattern counters.
std::function<std::size_t(const std::string&)> OccurrenceCounter(Storage t_storage, dret::Config& t_config);

// gcda.cpp, gcda_backend.cpp, dgcda.cpp, slp_ns.cpp
Built GCDA(Storage t_storage,
           dret::Config& t_config,
           uint32_t t_block_size,
           float t_storing_factor,
           axes::GCDASLPVariant t_slp,
           axes::GCDASetsCodec t_sets);
Built GCDABackend(Storage t_storage,
                  dret::Config& t_config,
                  uint32_t t_block_size,
                  float t_storing_factor,
                  axes::GetDocEnum t_get_doc,
                  axes::GCDASetsCodec t_sets);
Built DGCDA(Storage t_storage,
            dret::Config& t_config,
            uint32_t t_block_size,
            float t_storing_factor,
            axes::DGCDASLPVariant t_slp,
            axes::GCDASetsCodec t_sets);
Built SLPNS(Storage t_storage,
            dret::Config& t_config,
            axes::BareSLPVariant t_slp,
            uint32_t t_spacing,
            std::optional<axes::GCDASLPVariant> t_sampled,
            float t_storing_factor);

// rmq_<core>.cpp
struct RMQArgs {
  uint32_t block_size;
  float storing_factor;
  axes::GetDocEnum get_doc;
  axes::GCDASLPVariant gcda_slp;
  axes::BareSLPVariant bare_slp;
  axes::DGCDASLPVariant dgcda_slp;
  axes::RunValuesVariant run_values;  // ILCP-V, CILCP-V
  axes::PrevDocVariant prev_doc;      // SADA-V
  std::size_t sa_sampling;
};
Built Sada(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args);
Built Ilcp(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args);
Built Cilcp(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args);
Built SadaValue(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args);
Built IlcpValue(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args);
Built CilcpValue(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args);

// pdl_<codec>.cpp
struct PDLArgs {
  uint32_t block_size;
  float storing_factor;
  axes::GetDocEnum get_doc;
  axes::PDLStoragePolicy policy;
  axes::GCDASLPVariant gcda_slp;
  axes::BareSLPVariant bare_slp;
  axes::DGCDASLPVariant dgcda_slp;
  std::size_t sa_sampling;
  uint32_t backend_block_size;
};
Built PDLPlain(Storage t_storage, dret::Config& t_config, const PDLArgs& t_args);
Built PDLRP(Storage t_storage, dret::Config& t_config, const PDLArgs& t_args);
Built PDLBC(Storage t_storage, dret::Config& t_config, const PDLArgs& t_args);

}  // namespace bench::factories::query
