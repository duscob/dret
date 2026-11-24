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
#include <stdexcept>
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
#include "dret/backend/da.h"
#include "dret/backend/grammar.h"
#include "dret/backend/rlcsa.h"
#include "dret/backend/sa_phi.h"

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
using slp_ns::BareSLP_Diff;
using slp_ns::BareSLP_DiffEV;
using slp_ns::BareSLP_DiffDV;
using slp_ns::BareSLP_DiffVV;

// GetDoc policy template aliases. The GetDocSLP TSLP default intentionally
// matches gcda::SLP_Light, and GetDocDSLP's TSLP default matches the DGCDA
// default — cache sharing depends on the exact template match.
template <typename TStorage>
using GetDocDA = dret::backend::PlainDA<TStorage, kWidth>;

template <typename TStorage, typename TSLP = SLP_Light>
using GetDocSLP = dret::backend::GCDAGrammar<TStorage, kWidth, TSLP>;

template <typename TStorage, typename TSLP = BareSLP_IV>
using GetDocSLP_NS = dret::backend::BareGrammar<TStorage, kWidth, TSLP>;

// TDSLP defaults to the DGCDA default so the kDSLP cache is shared with the
// DGCDA build; other DGCDA variants (OTF/CRL/EV/DV/VV) select their own cache.
template <typename TStorage, typename TDSLP = dgcda::SLP_Default>
using GetDocDSLP = dret::backend::DiffGCDAGrammar<TStorage, kWidth, TDSLP>;

// SA-Phi backings — RLCSA-style baseline. The dret::backend::SAPhi class is
// templated on the locate-index type so the same implementation serves both
// r-index (dense, no sampling) and sr-index (subsampled, runtime sa_sampling).
template <typename TStorage>
using GetDocSAPhi_R  = dret::backend::SAPhi<
    TStorage, kWidth, sri::RIndex<TStorage, dret::Alphabet<>>>;
template <typename TStorage>
using GetDocSAPhi_SR = dret::backend::SAPhi<
    TStorage, kWidth, sri::SrIndexValidArea<TStorage, dret::Alphabet<>>>;

// RLCSA backing — paper-faithful batched locate(range) + getSequenceForPosition.
template <typename TStorage>
using GetDocRLCSA = dret::backend::RLCSA<TStorage, kWidth>;

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

// The value-rule cores (labelled -V in benchmark output, see axes.h): same runs
// as the published cores above, plus the array their value-based recursion stop
// consults. SADA-V persists prev_doc; ILCP-V and CILCP-V persist run_values.
// CILCP-V does NOT have its own partition: it and CilcpCore both build the CMR20
// runs through one shared routine and their run_heads are byte-identical,
// which CilcpPartitionTest asserts. They differ only in whether the per-run
// minima are stored.

// TPrevDoc container choices for the value-rule SADA-V core. Independent of
// TGetDoc; controls
// how the persisted per-SA-position prev_doc array is encoded on disk.
using PrevDoc_IV = sdsl::int_vector<>;   // fixed-width, bit-compressed (default)
using PrevDoc_DV = sdsl::dac_vector<>;
using PrevDoc_VV = sdsl::vlc_vector<>;

template <typename TStorage,
          typename TGetDoc = GetDocDA<TStorage>,
          typename TPrevDoc = PrevDoc_IV>
using SadaValueCore = dret::rmq::SadaValueCore<TStorage,
                                        kWidth,
                                        sdsl::rmq_succinct_sct<true>,
                                        sdsl::sd_vector<>,
                                        TGetDoc,
                                        TPrevDoc>;

// TRunValues container choices for the value-rule cores. Independent of TGetDoc;
// controls how the persisted per-run min(VILCP) array is encoded on disk.
using RunValues_IV = sdsl::int_vector<>;   // fixed-width, bit-compressed
using RunValues_DV = sdsl::dac_vector<>;   // direct access codes (default)
using RunValues_VV = sdsl::vlc_vector<>;   // variable-length codes

template <typename TStorage,
          typename TGetDoc = GetDocDA<TStorage>,
          typename TRunValues = RunValues_DV>
using IlcpValueCore = dret::rmq::IlcpValueCore<TStorage,
                                        kWidth,
                                        sdsl::sd_vector<>,
                                        sdsl::rmq_succinct_sct<true>,
                                        sdsl::sd_vector<>,
                                        TGetDoc,
                                        TRunValues>;

template <typename TStorage,
          typename TGetDoc = GetDocDA<TStorage>,
          typename TRunValues = RunValues_DV>
using CilcpValueCore = dret::rmq::CilcpValueCore<TStorage,
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

// Which core to build (2026-11 names, spec_version 2).
// SADA / ILCP / CILCP are the published algorithms (Sadakane 2007; Gagie et
// al. 2017; Cobas, Makinen, Rossi SPIRE 2020): they stop the recursion by
// marking reported documents and store no values. The -V cores are dret's
// value-rule variants: they also store the array the RMQ was built over and
// stop on its values.
enum class CoreKind { SADA_V, ILCP_V, CILCP_V, SADA, ILCP, CILCP };

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

// Build one RMQ-SLP-NS variant for a specific bare-SLP choice. No storing
// factor; the block size is the sample spacing of a differential TSLP and is
// ignored by a plain one (see GetDocSLP_NS).
template <template <typename, typename> class TCoreT, typename TStorage, typename TSLP>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSLP_NS(TStorage t_storage, dret::Config& t_config, uint32_t t_block_size) {
  using TCore = TCoreT<TStorage, GetDocSLP_NS<TStorage, TSLP>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, 0.0f);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename> class TCoreT, typename TStorage,
          typename TDSLP = dgcda::SLP_Default>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeDSLP(TStorage t_storage, dret::Config& t_config,
         uint32_t t_block_size, float t_storing_factor) {
  using TCore = TCoreT<TStorage, GetDocDSLP<TStorage, TDSLP>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, t_storing_factor);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

// Build one RMQ SA-Phi-R variant (r-index-backed; no sampling parameter).
// The published cores' constructors accept either 1 arg (storage) or 3 args
// (storage, bs, sf); for SA-Phi the inner GetDocSAPhi needs no per-core
// (bs, sf), so we use the 1-arg form — the inner GetDoc is
// default-constructed and its sa_sampling stays at 0 (the r-index sentinel,
// meaning "no subsampling" — i.e. dense sampling).
template <template <typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSAPhi_R(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocSAPhi_R<TStorage>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

// Build one RMQ-RLCSA variant (paper-faithful batched locate; no sampling knob).
// Same shape as MakeSAPhi_R — 1-arg constructor (storage only).
template <template <typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeRLCSA(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocRLCSA<TStorage>>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

// SA-Phi-SR (sr-index, subsampled) is deferred — see
// docs/pdl_rlcsa_baseline_plan.md. The user chose sa_phi_r-only for this
// stage. When SR is wired in, it will need either core-constructor
// extensions to accept the sampling rate, or a separate
// default-construct-then-set-sampling pattern on the inner GetDoc.

template <template <typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOne(TStorage t_storage, dret::Config& t_config,
        uint32_t t_block_size, float t_storing_factor,
        bench::axes::GetDocEnum t_get_doc,
        bench::axes::GCDASLPVariant t_gcda_slp,
        bench::axes::BareSLPVariant t_bare_slp,
        bench::axes::DGCDASLPVariant t_dgcda_slp,
        std::size_t t_sa_sampling = 0) {
  using bench::axes::GetDocEnum;
  using bench::axes::GCDASLPVariant;
  using bench::axes::BareSLPVariant;
  using bench::axes::DGCDASLPVariant;
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
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_Raw>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_DV>(t_storage, t_config, t_block_size);
        case BareSLPVariant::VV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_VV>(t_storage, t_config, t_block_size);
        case BareSLPVariant::Diff:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_Diff>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DiffEV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_DiffEV>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DiffDV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_DiffDV>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DiffVV:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_DiffVV>(t_storage, t_config, t_block_size);
        case BareSLPVariant::IV:
        default:
          return MakeSLP_NS<TCoreT, TStorage, BareSLP_IV>(t_storage, t_config, t_block_size);
      }
    case GetDocEnum::DSLP:
      switch (t_dgcda_slp) {
        case DGCDASLPVariant::OTF: return MakeDSLP<TCoreT, TStorage, dgcda::SLP_OTF>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::CRL: return MakeDSLP<TCoreT, TStorage, dgcda::SLP_CRL>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::EV:  return MakeDSLP<TCoreT, TStorage, dgcda::SLP_EV>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::DV:  return MakeDSLP<TCoreT, TStorage, dgcda::SLP_DV>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::VV:  return MakeDSLP<TCoreT, TStorage, dgcda::SLP_VV>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::Default:
        default:                   return MakeDSLP<TCoreT, TStorage, dgcda::SLP_Default>(t_storage, t_config, t_block_size, t_storing_factor);
      }
    case GetDocEnum::SAPhiR:
      return MakeSAPhi_R<TCoreT, TStorage>(t_storage, t_config);
    case GetDocEnum::SAPhiSR:
      // SA-Phi-SR deferred — fall through to SAPhi-R if requested.
      return MakeSAPhi_R<TCoreT, TStorage>(t_storage, t_config);
    case GetDocEnum::RLCSA: {
      // RMQ × RLCSA is INCOMPATIBLE: RLCSA builds a GSA-style SA where each
      // \x00 sequence separator gets a unique sort value (see utils.cpp
      // simpleSuffixSort line 385), while dret's SA treats all \x01 separators
      // as identical and compares across them. The two SAs cover the same
      // user-suffix SET but in different orders, so dret's RMQ argmin doesn't
      // map to the same suffix as RLCSA's compact[k-1-D]. SADA/ILCP miss docs
      // when the marker-based early stop triggers on the swapped doc. CILCP
      // happens to work (no early stop) but inefficiently. See
      // docs/pdl_rlcsa_baseline_plan.md for details. Fall through to DA.
      [[maybe_unused]] static const bool warned = [] {
        std::cerr << "[factories/rmq] RMQ × RLCSA is unsupported (SA-ordering "
                     "mismatch); falling through to DA. Use PDL × RLCSA for the "
                     "paper-faithful baseline.\n";
        return true;
      }();
      return MakeDA<TCoreT, TStorage>(t_storage, t_config);
    }
    case GetDocEnum::DA:
    default:
      return MakeDA<TCoreT, TStorage>(t_storage, t_config);
  }
}

//~~~~~~~  Value-rule-core dispatch (3-arg TCoreT taking <TStorage, TGetDoc, TRunValues>)
//
// IlcpValueCore / CilcpValueCore take an extra TRunValues template parameter on top
// of the standard 2-arg core template, so the existing 2-arg MakeOne path
// doesn't fit. The MakeXxx_S helpers mirror MakeXxx but plumb TRunValues
// through; MakeOnePub_T dispatches on TGetDoc with TRunValues fixed; MakeOneValue
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

// SA-Phi-R needs no sampling parameter, so the storage-only constructor builds
// it as for the published cores (see MakeSAPhi_R).
template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeSAPhi_R_S(TStorage t_storage, dret::Config& t_config) {
  using TCore = TCoreT<TStorage, GetDocSAPhi_R<TStorage>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage);
  auto idx = std::make_shared<TIndex>(t_storage, core);
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
MakeSLP_NS_S(TStorage t_storage, dret::Config& t_config, uint32_t t_block_size) {
  using TCore = TCoreT<TStorage, GetDocSLP_NS<TStorage, TSLP>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, 0.0f);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

template <template <typename, typename, typename> class TCoreT, typename TStorage,
          typename TRunValues, typename TDSLP = dgcda::SLP_Default>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeDSLP_S(TStorage t_storage, dret::Config& t_config,
           uint32_t t_block_size, float t_storing_factor) {
  using TCore = TCoreT<TStorage, GetDocDSLP<TStorage, TDSLP>, TRunValues>;
  using TIndex = Idx<TStorage, TCore>;
  TCore core(t_storage, t_block_size, t_storing_factor);
  auto idx = std::make_shared<TIndex>(t_storage, core);
  construct(*idx, t_config);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

// NOTE: the value-rule cores (IlcpLikeValueCore / SadaValueCore) support SA-Phi-R
// only: their 1-arg constructors cannot pass the sa_sampling rate SA-Phi-SR
// needs, so SA-Phi-SR builds SA-Phi-R, as it does for the published cores. RLCSA
// is unsupported and throws rather than measuring another backend under its
// name (it used to fall through to the DA).

template <template <typename, typename, typename> class TCoreT, typename TStorage, typename TRunValues>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOnePub_T(TStorage t_storage, dret::Config& t_config,
           uint32_t t_block_size, float t_storing_factor,
           bench::axes::GetDocEnum t_get_doc,
           bench::axes::GCDASLPVariant t_gcda_slp,
           bench::axes::BareSLPVariant t_bare_slp,
           bench::axes::DGCDASLPVariant t_dgcda_slp,
           std::size_t t_sa_sampling = 0) {
  using bench::axes::GetDocEnum;
  using bench::axes::GCDASLPVariant;
  using bench::axes::BareSLPVariant;
  using bench::axes::DGCDASLPVariant;
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
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_Raw, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_DV, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::VV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_VV, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::Diff:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_Diff, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DiffEV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_DiffEV, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DiffDV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_DiffDV, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::DiffVV:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_DiffVV, TRunValues>(t_storage, t_config, t_block_size);
        case BareSLPVariant::IV:
        default:
          return MakeSLP_NS_S<TCoreT, TStorage, BareSLP_IV, TRunValues>(t_storage, t_config, t_block_size);
      }
    case GetDocEnum::DSLP:
      switch (t_dgcda_slp) {
        case DGCDASLPVariant::OTF: return MakeDSLP_S<TCoreT, TStorage, TRunValues, dgcda::SLP_OTF>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::CRL: return MakeDSLP_S<TCoreT, TStorage, TRunValues, dgcda::SLP_CRL>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::EV:  return MakeDSLP_S<TCoreT, TStorage, TRunValues, dgcda::SLP_EV>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::DV:  return MakeDSLP_S<TCoreT, TStorage, TRunValues, dgcda::SLP_DV>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::VV:  return MakeDSLP_S<TCoreT, TStorage, TRunValues, dgcda::SLP_VV>(t_storage, t_config, t_block_size, t_storing_factor);
        case DGCDASLPVariant::Default:
        default:                   return MakeDSLP_S<TCoreT, TStorage, TRunValues, dgcda::SLP_Default>(t_storage, t_config, t_block_size, t_storing_factor);
      }
    case GetDocEnum::SAPhiR:
    case GetDocEnum::SAPhiSR:
      return MakeSAPhi_R_S<TCoreT, TStorage, TRunValues>(t_storage, t_config);
    case GetDocEnum::RLCSA:
      throw std::invalid_argument("RMQ: the value-rule cores do not support get_doc=rlcsa");
    case GetDocEnum::DA:
    default:
      return MakeDA_S<TCoreT, TStorage, TRunValues>(t_storage, t_config);
  }
}

template <template <typename, typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOneValue(TStorage t_storage, dret::Config& t_config,
         uint32_t t_block_size, float t_storing_factor,
         bench::axes::GetDocEnum t_get_doc,
         bench::axes::GCDASLPVariant t_gcda_slp,
         bench::axes::BareSLPVariant t_bare_slp,
         bench::axes::DGCDASLPVariant t_dgcda_slp,
         bench::axes::RunValuesVariant t_run_values,
         std::size_t t_sa_sampling = 0) {
  using bench::axes::RunValuesVariant;
  switch (t_run_values) {
    case RunValuesVariant::IV:
      return MakeOnePub_T<TCoreT, TStorage, RunValues_IV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
    case RunValuesVariant::VV:
      return MakeOnePub_T<TCoreT, TStorage, RunValues_VV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
    case RunValuesVariant::DV:
    default:
      return MakeOnePub_T<TCoreT, TStorage, RunValues_DV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
  }
}

// The published SADA core has TPrevDoc instead of TRunValues. Mirrors MakeOneValue but
// dispatches on PrevDocVariant; the inner T-dispatch reuses the existing
// MakeXxx_S helpers (TCoreT here is a 3-arg <TStorage, TGetDoc, TPrevDoc>
// template — same shape as IlcpValueCore / CilcpValueCore).
template <template <typename, typename, typename> class TCoreT, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeOneSadaValue(TStorage t_storage, dret::Config& t_config,
              uint32_t t_block_size, float t_storing_factor,
              bench::axes::GetDocEnum t_get_doc,
              bench::axes::GCDASLPVariant t_gcda_slp,
              bench::axes::BareSLPVariant t_bare_slp,
              bench::axes::DGCDASLPVariant t_dgcda_slp,
              bench::axes::PrevDocVariant t_prev_doc,
              std::size_t t_sa_sampling = 0) {
  using bench::axes::PrevDocVariant;
  switch (t_prev_doc) {
    case PrevDocVariant::DV:
      return MakeOnePub_T<TCoreT, TStorage, PrevDoc_DV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
    case PrevDocVariant::VV:
      return MakeOnePub_T<TCoreT, TStorage, PrevDoc_VV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
    case PrevDocVariant::IV:
    default:
      return MakeOnePub_T<TCoreT, TStorage, PrevDoc_IV>(
          t_storage, t_config, t_block_size, t_storing_factor, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
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
     bench::axes::DGCDASLPVariant t_dgcda_slp = bench::axes::DGCDASLPVariant::Default,
     bench::axes::RunValuesVariant t_run_values = bench::axes::RunValuesVariant::DV,
     bench::axes::PrevDocVariant t_prev_doc = bench::axes::PrevDocVariant::IV,
     std::size_t t_sa_sampling = 0) {
  switch (t_core) {
    case CoreKind::ILCP:
      return detail::MakeOne<IlcpCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                        t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
    case CoreKind::CILCP:
      return detail::MakeOne<CilcpCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                         t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
    case CoreKind::SADA_V:
      return detail::MakeOneSadaValue<SadaValueCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                               t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_prev_doc, t_sa_sampling);
    case CoreKind::ILCP_V:
      return detail::MakeOneValue<IlcpValueCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                          t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_run_values, t_sa_sampling);
    case CoreKind::CILCP_V:
      return detail::MakeOneValue<CilcpValueCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                           t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_run_values, t_sa_sampling);
    case CoreKind::SADA:
    default:
      return detail::MakeOne<SadaCore>(t_storage, t_config, t_block_size, t_storing_factor,
                                        t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp, t_sa_sampling);
  }
}

}  // namespace bench::factories::rmq
