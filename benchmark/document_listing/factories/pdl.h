//
// PDL factory header: the types of the dret::pdl::DocListIdxPDL family and
// ByGetDoc, which makes one cell's index. pdl.cpp turns the sweep into cells;
// each codec's are made in pdl_<codec>.cpp.
//
// The PDL family is a 2D type-level axis (codec × get-doc backing) times a
// runtime storage-policy axis. Codec choice picks one of the three
// DocListIdxPDL{Plain,RP,BC} alias templates; get-doc choice picks the
// raw-range backing (DA / SLP / DSLP). Storage policy is threaded into the
// constructor as a runtime argument.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <sdsl/util.hpp>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/backend/da.h"
#include "dret/backend/grammar.h"
#include "dret/backend/rlcsa.h"
#include "dret/backend/sa_phi.h"
#include "dret/pdl/set_codecs.h"
#include "dret/pdl/storage_policy.h"

#include "../axes.h"
#include "../enum_traits.h"
#include "cells.h"
#include "dgcda.h"
#include "gcda.h"
#include "pdl_codecs.h"
#include "slp_ns.h"

namespace bench::factories::pdl {

inline constexpr std::uint8_t kWidth = dret::Alphabet<>::int_width;

// GetDoc backing aliases. SLP / SLP-NS / DSLP are parameterised on the grammar
// type so PDL can pick any GCDA-family backend (the same TSLP defaults the RMQ
// and GCDA/DGCDA/SLP-NS factories use, so the typed cache files are shared).
template <typename TStorage>
using GetDocsDA   = dret::backend::PlainDA<TStorage, kWidth>;

template <typename TStorage, typename TSLP = gcda::SLP_Light>
using GetDocsSLP  = dret::backend::GCDAGrammar<TStorage, kWidth, TSLP>;

template <typename TStorage, typename TSLP = slp_ns::BareSLP_IV>
using GetDocsSLP_NS = dret::backend::BareGrammar<TStorage, kWidth, TSLP>;

template <typename TStorage, typename TDSLP = dgcda::SLP_Default>
using GetDocsDSLP = dret::backend::DiffGCDAGrammar<TStorage, kWidth, TDSLP>;

// SA-Phi backings — the RLCSA-style PDL baseline. R variant uses dense
// r-index samples (no sampling knob); SR variant uses subsampled sr-index
// samples (sa_sampling knob threaded through the spec).
template <typename TStorage>
using GetDocsSAPhi_R  = dret::backend::SAPhi<TStorage, kWidth>;

template <typename TStorage>
using GetDocsSAPhi_SR = dret::backend::SAPhi<TStorage, kWidth, sri::SrIndexValidArea<TStorage, dret::Alphabet<kWidth>>>;

// RLCSA backing — paper-faithful batched locate(range) + getSequenceForPosition.
template <typename TStorage>
using GetDocsRLCSA = dret::backend::RLCSA<TStorage, kWidth>;

// Typed-index aliases. One template per codec; TGetDocs defaults to DA.
// TStorage and TGetDocs are both parametric.
template <typename TStorage, typename TGetDocs = GetDocsDA<TStorage>>
using IdxPlain = dret::pdl::DocListIdxPDLPlain<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TGetDocs>;

template <typename TStorage, typename TGetDocs = GetDocsDA<TStorage>>
using IdxRP = dret::pdl::DocListIdxPDLRP<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TGetDocs>;

template <typename TStorage, typename TGetDocs = GetDocsDA<TStorage>>
using IdxBC = dret::pdl::DocListIdxPDLBC<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TGetDocs>;

namespace detail {

// The index over its default backend, which it builds from PDL's own block size.
template <typename TIndex>
MakeIndex OverDefaultBackend(const Params& t_p) {
  return MakeOf<TIndex>(t_p.block_size, t_p.storing_factor, axes::toPDLStoragePolicy(t_p.policy));
}

// The index over a backend with its own block size, t_p.backend_block_size, independent of PDL's
// (block size, storing factor), as RMQ builds the same backend from its block size (so both families
// share its cache files). The default constructor would build the backend from PDL's block size.
template <typename TIndex, typename TGetDocs>
MakeIndex OverBackend(const Params& t_p, float t_backend_storing_factor) {
  return MakeBy<TIndex>([t_p, t_backend_storing_factor](Storage t_storage) {
    TGetDocs backend(t_storage, t_p.backend_block_size, t_backend_storing_factor);
    return std::make_shared<TIndex>(
        t_storage, backend, t_p.block_size, t_p.storing_factor, axes::toPDLStoragePolicy(t_p.policy));
  });
}

}  // namespace detail

// The index of one PDL cell with the stored-set codec TCodec, over the get-doc backend the cell
// names (and its grammar variant). GCDA-light and the differential bare grammars are built apart, at
// the backend block size: GCDA-light with storing factor 32 (it only shapes GCDA's document lists,
// which the backend does not load), a differential grammar sampled every backend-block-size
// positions. SA-Phi-SR builds SA-Phi-R (the sr-index backend is not wired in).
template <template <typename, typename> class TCodec>
MakeIndex ByGetDoc(const Params& t_p) {
  using namespace bench::axes;
  using detail::OverBackend;
  using detail::OverDefaultBackend;
  switch (t_p.get_doc) {
    case GetDocEnum::SLP:
      switch (t_p.gcda_slp) {
        case GCDASLPVariant::CompactBP:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP<Storage, gcda::SLP_CompactBP>>>(t_p);
        case GCDASLPVariant::CompactLOUDS:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP<Storage, gcda::SLP_CompactLOUDS>>>(t_p);
        case GCDASLPVariant::Combined:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP<Storage, gcda::SLP_Combined>>>(t_p);
        case GCDASLPVariant::Light:
        default: {
          using TGD = GetDocsSLP<Storage, gcda::SLP_Light>;
          return OverBackend<TCodec<Storage, TGD>, TGD>(t_p, 32.0f);
        }
      }
    case GetDocEnum::SLP_NS:
      switch (t_p.bare_slp) {
        case BareSLPVariant::Raw:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP_NS<Storage, slp_ns::BareSLP_Raw>>>(t_p);
        case BareSLPVariant::DV:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP_NS<Storage, slp_ns::BareSLP_DV>>>(t_p);
        case BareSLPVariant::VV:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP_NS<Storage, slp_ns::BareSLP_VV>>>(t_p);
        case BareSLPVariant::Diff: {
          using TGD = GetDocsSLP_NS<Storage, slp_ns::BareSLP_Diff>;
          return OverBackend<TCodec<Storage, TGD>, TGD>(t_p, 0.0f);
        }
        case BareSLPVariant::DiffEV: {
          using TGD = GetDocsSLP_NS<Storage, slp_ns::BareSLP_DiffEV>;
          return OverBackend<TCodec<Storage, TGD>, TGD>(t_p, 0.0f);
        }
        case BareSLPVariant::DiffDV: {
          using TGD = GetDocsSLP_NS<Storage, slp_ns::BareSLP_DiffDV>;
          return OverBackend<TCodec<Storage, TGD>, TGD>(t_p, 0.0f);
        }
        case BareSLPVariant::DiffVV: {
          using TGD = GetDocsSLP_NS<Storage, slp_ns::BareSLP_DiffVV>;
          return OverBackend<TCodec<Storage, TGD>, TGD>(t_p, 0.0f);
        }
        case BareSLPVariant::IV:
        default:
          return OverDefaultBackend<TCodec<Storage, GetDocsSLP_NS<Storage, slp_ns::BareSLP_IV>>>(t_p);
      }
    case GetDocEnum::DSLP:
      switch (t_p.dgcda_slp) {
        case DGCDASLPVariant::OTF:
          return OverDefaultBackend<TCodec<Storage, GetDocsDSLP<Storage, dgcda::SLP_OTF>>>(t_p);
        case DGCDASLPVariant::CRL:
          return OverDefaultBackend<TCodec<Storage, GetDocsDSLP<Storage, dgcda::SLP_CRL>>>(t_p);
        case DGCDASLPVariant::EV:
          return OverDefaultBackend<TCodec<Storage, GetDocsDSLP<Storage, dgcda::SLP_EV>>>(t_p);
        case DGCDASLPVariant::DV:
          return OverDefaultBackend<TCodec<Storage, GetDocsDSLP<Storage, dgcda::SLP_DV>>>(t_p);
        case DGCDASLPVariant::VV:
          return OverDefaultBackend<TCodec<Storage, GetDocsDSLP<Storage, dgcda::SLP_VV>>>(t_p);
        case DGCDASLPVariant::Default:
        default:
          return OverDefaultBackend<TCodec<Storage, GetDocsDSLP<Storage, dgcda::SLP_Default>>>(t_p);
      }
    case GetDocEnum::SAPhiR:
    case GetDocEnum::SAPhiSR:
      return OverDefaultBackend<TCodec<Storage, GetDocsSAPhi_R<Storage>>>(t_p);
    case GetDocEnum::RLCSA:
      return OverDefaultBackend<TCodec<Storage, GetDocsRLCSA<Storage>>>(t_p);
    case GetDocEnum::DA:
    default:
      return OverDefaultBackend<TCodec<Storage, GetDocsDA<Storage>>>(t_p);
  }
}

}  // namespace bench::factories::pdl
