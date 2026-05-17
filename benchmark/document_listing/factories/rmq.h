//
// RMQ factory header — type aliases + Make() for the SADA / ILCP / CILCP
// document-listing cores in dret::rmq.
//
// The RMQ family has the largest cross-product of any doc-list family:
//   3 cores  ×  {DA, SLP, SLP-NS, DSLP}  ×  (per-get-doc TSLP sub-axis).
// We expose four GetDoc template aliases, three Core template aliases (each
// parameterised on TStorage and TGetDoc), and one Idx template alias. This
// keeps the type-alias count bounded and lets callers compose the exact
// instantiation they want.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <sdsl/dac_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/vlc_vector.hpp>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/rmq/rmq_get_doc_policies.h"

#include "../axes.h"
#include "dgcda.h"
#include "gcda.h"
#include "slp_ns.h"

namespace bench::factories::rmq {

inline constexpr std::uint8_t kWidth = dret::Alphabet<>::int_width;

// Re-export the sibling-factory SLP type aliases so RMQ-SLP / RMQ-SLP-NS /
// RMQ-DSLP cache files match the corresponding GCDA / SLP-NS / DGCDA builds
// via grammar:: type-hashing on disk.
using gcda::SLP_Light;
using gcda::SLP_CompactBP;
using gcda::SLP_CompactLOUDS;
using gcda::SLP_Combined;
using slp_ns::BareSLP_IV;
using slp_ns::BareSLP_Raw;
using slp_ns::BareSLP_DV;
using slp_ns::BareSLP_VV;

// GetDoc policy template aliases. The GetDocSLP TSLP default intentionally
// matches gcda::SLP_Light, and GetDocDSLP's TSLP default matches the DGCDA
// default — cache sharing depends on the exact template match.
template <typename TStorage>
using GetDocDA = dret::rmq::GetDocDA<TStorage, kWidth>;

template <typename TStorage, typename TSLP = SLP_Light>
using GetDocSLP = dret::rmq::GetDocSLP<TStorage, kWidth, TSLP>;

template <typename TStorage, typename TSLP = BareSLP_IV>
using GetDocSLP_NS = dret::rmq::GetDocSLP_NS<TStorage, kWidth, TSLP>;

template <typename TStorage>
using GetDocDSLP = dret::rmq::GetDocDSLP<TStorage, kWidth>;

// Core template aliases. TGetDoc defaults to the DA backing.
template <typename TStorage, typename TGetDoc = GetDocDA<TStorage>>
using SadaCore = dret::rmq::SadaCore<TStorage,
                                      kWidth,
                                      sdsl::rmq_succinct_sct<true>,
                                      sdsl::sd_vector<>,
                                      TGetDoc>;

template <typename TStorage, typename TGetDoc = GetDocDA<TStorage>>
using IlcpCore = dret::rmq::IlcpCore<TStorage,
                                      kWidth,
                                      sdsl::sd_vector<>,
                                      sdsl::rmq_succinct_sct<true>,
                                      sdsl::sd_vector<>,
                                      TGetDoc>;

template <typename TStorage, typename TGetDoc = GetDocDA<TStorage>>
using CilcpCore = dret::rmq::CilcpCore<TStorage,
                                        kWidth,
                                        sdsl::sd_vector<>,
                                        sdsl::rmq_succinct_sct<true>,
                                        sdsl::sd_vector<>,
                                        TGetDoc>;

// Sadakane-style (-S) parallel cores: same RMQ data, canonical depth-based
// recursion stop. SADA-S persists prev_doc; ILCP-S and CILCP-S persist
// run_values. CILCP-S's RLE follows paper Def 1 (CILCP★) and is distinct
// from existing CilcpCore.
template <typename TStorage, typename TGetDoc = GetDocDA<TStorage>>
using SadaSCore = dret::rmq::SadaSCore<TStorage,
                                        kWidth,
                                        sdsl::rmq_succinct_sct<true>,
                                        sdsl::sd_vector<>,
                                        TGetDoc>;

// TRunValues container choices for the -S families. Independent of TGetDoc;
// controls how the persisted per-run min(VILCP) array is encoded on disk.
using RunValues_IV = sdsl::int_vector<>;   // fixed-width, bit-compressed
using RunValues_DV = sdsl::dac_vector<>;   // direct access codes (default)
using RunValues_VV = sdsl::vlc_vector<>;   // variable-length codes

template <typename TStorage,
          typename TGetDoc = GetDocDA<TStorage>,
          typename TRunValues = RunValues_DV>
using IlcpSCore = dret::rmq::IlcpSCore<TStorage,
                                        kWidth,
                                        sdsl::sd_vector<>,
                                        sdsl::rmq_succinct_sct<true>,
                                        sdsl::sd_vector<>,
                                        TGetDoc,
                                        TRunValues>;

template <typename TStorage,
          typename TGetDoc = GetDocDA<TStorage>,
          typename TRunValues = RunValues_DV>
using CilcpSCore = dret::rmq::CilcpSCore<TStorage,
                                          kWidth,
                                          sdsl::sd_vector<>,
                                          sdsl::rmq_succinct_sct<true>,
                                          sdsl::sd_vector<>,
                                          TGetDoc,
                                          TRunValues>;

// Typed-index template alias. Storage and Core are both parametric.
template <typename TStorage, typename TCore>
using Idx = dret::rmq::DocListIdxRMQ<TStorage,
                                      dret::Alphabet<>,
                                      sri::RIndexCount<TStorage, dret::Alphabet<>>,
                                      TCore>;

// Which core to build.
enum class CoreKind { SADA, ILCP, CILCP, SADA_S, ILCP_S, CILCP_S };

namespace detail {

// Build one RMQ-DA variant (no bs/sf knobs).
template <template <typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeDA(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocDA<TStorage>>;
  using TIndex = Idx<TStorage, TCore>;
  auto idx = std::make_shared<TIndex>(t_storage);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

// Build one RMQ-SLP variant for a specific GCDA TSLP choice. (bs, sf) thread
// into the GetDoc-SLP cache so the SLP file is shared with the corresponding
// GCDA build via grammar:: type-hashing on disk.
template <template <typename, typename> class TCoreT, typename TStorage, typename TSLP>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSLP(TStorage t_storage, dret::Config& t_config,
        uint32_t t_block_size, float t_storing_factor) {
  using TCore = TCoreT<TStorage, GetDocSLP<TStorage, TSLP>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, t_storing_factor);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

// Build one RMQ-SLP-NS variant for a specific bare-SLP choice (no bs/sf).
template <template <typename, typename> class TCoreT, typename TStorage, typename TSLP>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSLP_NS(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocSLP_NS<TStorage, TSLP>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeDSLP(TStorage t_storage, dret::Config& t_config,
         uint32_t t_block_size, float t_storing_factor) {
  using TCore = TCoreT<TStorage, GetDocDSLP<TStorage>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, t_storing_factor);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOne(TStorage t_storage, dret::Config& t_config,
        uint32_t t_block_size, float t_storing_factor,
        bench::axes::GetDocEnum t_get_doc,
        bench::axes::GCDASLPVariant t_gcda_slp,
        bench::axes::BareSLPVariant t_bare_slp) {
  using bench::axes::GetDocEnum;
  using bench::axes::GCDASLPVariant;
  using bench::axes::BareSLPVariant;
  switch (t_get_doc) {
    case GetDocEnum::SLP:
      switch (t_gcda_slp) {
        case GCDASLPVariant::CompactBP:
          return MakeSLP<TCoreT, TStorage, SLP_CompactBP>(t_storage, t_config, t_block_size, t_storing_factor);
        case GCDASLPVariant::CompactLOUDS:
          return MakeSLP<TCoreT, TStorage, SLP_CompactLOUDS>(t_storage, t_config, t_block_size, t_storing_factor);
        case GCDASLPVariant::Combined:
          return MakeSLP<TCoreT, TStorage, SLP_Combined>(t_storage, t_config, t_block_size, t_storing_factor);
        case GCDASLPVariant::Light:
        default:
          return MakeSLP<TCoreT, TStorage, SLP_Light>(t_storage, t_config, t_block_size, t_storing_factor);
      }
    case GetDocEnum::SLP_NS:
      switch (t_bare_slp) {
        case BareSLPVariant::Raw:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_Raw>(t_storage, t_config);
        case BareSLPVariant::DV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_DV>(t_storage, t_config);
        case BareSLPVariant::VV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_VV>(t_storage, t_config);
        case BareSLPVariant::IV:
        default:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_IV>(t_storage, t_config);
      }
    case GetDocEnum::DSLP:
      return MakeDSLP<TCoreT, TStorage>(t_storage, t_config, t_block_size, t_storing_factor);
    case GetDocEnum::DA:
    default:
      return MakeDA<TCoreT, TStorage>(t_storage, t_config);
  }
}

//~~~~~~~  -S family dispatch (3-arg TCoreT taking <TStorage, TGetDoc, TRunValues>)
//
// IlcpSCore / CilcpSCore take an extra TRunValues template parameter on top
// of the standard 2-arg core template, so the existing 2-arg MakeOne path
// doesn't fit. The MakeXxx_S helpers mirror MakeXxx but plumb TRunValues
// through; MakeOneS_T dispatches on TGetDoc with TRunValues fixed; MakeOneS
// fans out across TRunValues.

template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeDA_S(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocDA<TStorage>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  auto idx = std::make_shared<TIndex>(t_storage);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TSLP, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSLP_S(TStorage t_storage, dret::Config& t_config,
          uint32_t t_block_size, float t_storing_factor) {
  using TCore = TCoreT<TStorage, GetDocSLP<TStorage, TSLP>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, t_storing_factor);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TSLP, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSLP_NS_S(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocSLP_NS<TStorage, TSLP>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeDSLP_S(TStorage t_storage, dret::Config& t_config,
           uint32_t t_block_size, float t_storing_factor) {
  using TCore = TCoreT<TStorage, GetDocDSLP<TStorage>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, t_storing_factor);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOneS_T(TStorage t_storage, dret::Config& t_config,
           uint32_t t_block_size, float t_storing_factor,
           bench::axes::GetDocEnum t_get_doc,
           bench::axes::GCDASLPVariant t_gcda_slp,
           bench::axes::BareSLPVariant t_bare_slp) {
  using bench::axes::GetDocEnum;
  using bench::axes::GCDASLPVariant;
  using bench::axes::BareSLPVariant;
  switch (t_get_doc) {
    case GetDocEnum::SLP:
      switch (t_gcda_slp) {
        case GCDASLPVariant::CompactBP:
          return MakeSLP_S<TCoreT, TStorage, SLP_CompactBP, TRunValues>(t_storage, t_config, t_block_size, t_storing_factor);
        case GCDASLPVariant::CompactLOUDS:
          return MakeSLP_S<TCoreT, TStorage, SLP_CompactLOUDS, TRunValues>(t_storage, t_config, t_block_size, t_storing_factor);
        case GCDASLPVariant::Combined:
          return MakeSLP_S<TCoreT, TStorage, SLP_Combined, TRunValues>(t_storage, t_config, t_block_size, t_storing_factor);
        case GCDASLPVariant::Light:
        default:
          return MakeSLP_S<TCoreT, TStorage, SLP_Light, TRunValues>(t_storage, t_config, t_block_size, t_storing_factor);
      }
    case GetDocEnum::SLP_NS:
      switch (t_bare_slp) {
        case BareSLPVariant::Raw:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_Raw, TRunValues>(t_storage, t_config);
        case BareSLPVariant::DV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_DV, TRunValues>(t_storage, t_config);
        case BareSLPVariant::VV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_VV, TRunValues>(t_storage, t_config);
        case BareSLPVariant::IV:
        default:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_IV, TRunValues>(t_storage, t_config);
      }
    case GetDocEnum::DSLP:
      return MakeDSLP_S<TCoreT, TStorage, TRunValues>(t_storage, t_config, t_block_size, t_storing_factor);
    case GetDocEnum::DA:
    default:
      return MakeDA_S<TCoreT, TStorage, TRunValues>(t_storage, t_config);
  }
}

template <template <typename, typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOneS(TStorage t_storage, dret::Config& t_config,
         uint32_t t_block_size, float t_storing_factor,
         bench::axes::GetDocEnum t_get_doc,
         bench::axes::GCDASLPVariant t_gcda_slp,
         bench::axes::BareSLPVariant t_bare_slp,
         bench::axes::RunValuesVariant t_run_values) {
  using bench::axes::RunValuesVariant;
  switch (t_run_values) {
    case RunValuesVariant::IV:
      return MakeOneS_T<TCoreT, TStorage, RunValues_IV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp);
    case RunValuesVariant::VV:
      return MakeOneS_T<TCoreT, TStorage, RunValues_VV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp);
    case RunValuesVariant::DV:
    default:
      return MakeOneS_T<TCoreT, TStorage, RunValues_DV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp);
  }
}

}  // namespace detail

template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage, dret::Config& t_config,
     uint32_t t_block_size, float t_storing_factor,
     CoreKind t_core,
     bench::axes::GetDocEnum t_get_doc,
     bench::axes::GCDASLPVariant t_gcda_slp,
     bench::axes::BareSLPVariant t_bare_slp,
     bench::axes::RunValuesVariant t_run_values = bench::axes::RunValuesVariant::DV) {
  switch (t_core) {
    case CoreKind::ILCP:
      return detail::MakeOne<IlcpCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                        t_get_doc, t_gcda_slp, t_bare_slp);
    case CoreKind::CILCP:
      return detail::MakeOne<CilcpCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                         t_get_doc, t_gcda_slp, t_bare_slp);
    case CoreKind::SADA_S:
      return detail::MakeOne<SadaSCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                         t_get_doc, t_gcda_slp, t_bare_slp);
    case CoreKind::ILCP_S:
      return detail::MakeOneS<IlcpSCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                          t_get_doc, t_gcda_slp, t_bare_slp, t_run_values);
    case CoreKind::CILCP_S:
      return detail::MakeOneS<CilcpSCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                           t_get_doc, t_gcda_slp, t_bare_slp, t_run_values);
    case CoreKind::SADA:
    default:
      return detail::MakeOne<SadaCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                        t_get_doc, t_gcda_slp, t_bare_slp);
  }
}

}  // namespace bench::factories::rmq
