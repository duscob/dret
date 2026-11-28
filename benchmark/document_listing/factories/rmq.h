//
// RMQ factory header: the types of the dret::rmq cores (SADA / ILCP / CILCP and
// their value-rule -V variants) and ByGetDoc, which makes one cell's index.
// rmq.cpp turns the sweep into cells; each core's are made in rmq_<core>.cpp.
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
#include "cells.h"
#include "dgcda.h"
#include "gcda.h"
#include "rmq_cores.h"
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

namespace detail {

// The index over its default core, which builds its own backend.
template <typename TCore>
MakeIndex OverDefaultCore() {
  return MakeOf<Idx<Storage, TCore>>();
}

// The index over a core made as TCore(storage, t_args...).
template <typename TCore, typename... TArgs>
MakeIndex OverCore(TArgs... t_args) {
  using TIndex = Idx<Storage, TCore>;
  return MakeBy<TIndex>(
      [=](Storage t_storage) { return std::make_shared<TIndex>(t_storage, TCore(t_storage, t_args...)); });
}

}  // namespace detail

// The index of one RMQ cell, its core TCoreOf<TGetDoc> over the backend the cell names. A grammar
// backend is built at the cell's (block size, storing factor), the same types the GCDA / DGCDA /
// SLP-NS factories use, so the grammar cache files are shared with theirs.
//
// SA-Phi-SR builds SA-Phi-R: the cores' constructors cannot pass the subsample rate SA-Phi-SR needs.
// RLCSA is not supported: RLCSA builds a GSA-style suffix array, where each sequence separator sorts
// apart, while dret's treats all separators as equal, so the RMQ argmin does not map to the suffix
// RLCSA locates and the marking cores miss documents (docs/pdl_rlcsa_baseline_plan.md). PDL over the
// RLCSA is the paper-faithful baseline.
template <template <typename> class TCoreOf>
MakeIndex ByGetDoc(const Params& t_p) {
  using namespace bench::axes;
  const auto bs = t_p.block_size;
  const auto sf = t_p.storing_factor;
  switch (t_p.get_doc) {
    case GetDocEnum::SLP:
      switch (t_p.gcda_slp) {
        case GCDASLPVariant::CompactBP:
          return detail::OverCore<TCoreOf<GetDocSLP<Storage, SLP_CompactBP>>>(bs, sf);
        case GCDASLPVariant::CompactLOUDS:
          return detail::OverCore<TCoreOf<GetDocSLP<Storage, SLP_CompactLOUDS>>>(bs, sf);
        case GCDASLPVariant::Combined:
          return detail::OverCore<TCoreOf<GetDocSLP<Storage, SLP_Combined>>>(bs, sf);
        case GCDASLPVariant::Light:
        default:
          return detail::OverCore<TCoreOf<GetDocSLP<Storage, SLP_Light>>>(bs, sf);
      }
    case GetDocEnum::SLP_NS:
      // No storing factor; the block size is the sample spacing of a differential grammar and is
      // ignored by a plain one.
      switch (t_p.bare_slp) {
        case BareSLPVariant::Raw:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_Raw>>>(bs, 0.0f);
        case BareSLPVariant::DV:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_DV>>>(bs, 0.0f);
        case BareSLPVariant::VV:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_VV>>>(bs, 0.0f);
        case BareSLPVariant::Diff:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_Diff>>>(bs, 0.0f);
        case BareSLPVariant::DiffEV:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_DiffEV>>>(bs, 0.0f);
        case BareSLPVariant::DiffDV:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_DiffDV>>>(bs, 0.0f);
        case BareSLPVariant::DiffVV:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_DiffVV>>>(bs, 0.0f);
        case BareSLPVariant::IV:
        default:
          return detail::OverCore<TCoreOf<GetDocSLP_NS<Storage, BareSLP_IV>>>(bs, 0.0f);
      }
    case GetDocEnum::DSLP:
      switch (t_p.dgcda_slp) {
        case DGCDASLPVariant::OTF:
          return detail::OverCore<TCoreOf<GetDocDSLP<Storage, dgcda::SLP_OTF>>>(bs, sf);
        case DGCDASLPVariant::CRL:
          return detail::OverCore<TCoreOf<GetDocDSLP<Storage, dgcda::SLP_CRL>>>(bs, sf);
        case DGCDASLPVariant::EV:
          return detail::OverCore<TCoreOf<GetDocDSLP<Storage, dgcda::SLP_EV>>>(bs, sf);
        case DGCDASLPVariant::DV:
          return detail::OverCore<TCoreOf<GetDocDSLP<Storage, dgcda::SLP_DV>>>(bs, sf);
        case DGCDASLPVariant::VV:
          return detail::OverCore<TCoreOf<GetDocDSLP<Storage, dgcda::SLP_VV>>>(bs, sf);
        case DGCDASLPVariant::Default:
        default:
          return detail::OverCore<TCoreOf<GetDocDSLP<Storage, dgcda::SLP_Default>>>(bs, sf);
      }
    case GetDocEnum::SAPhiR:
    case GetDocEnum::SAPhiSR:
      return detail::OverCore<TCoreOf<GetDocSAPhi_R<Storage>>>();
    case GetDocEnum::RLCSA:
      throw std::invalid_argument("rmq: get_doc=rlcsa is not supported (see factories/rmq.h)");
    case GetDocEnum::DA:
    default:
      return detail::OverDefaultCore<TCoreOf<GetDocDA<Storage>>>();
  }
}

}  // namespace bench::factories::rmq
