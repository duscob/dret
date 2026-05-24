//
// PDL factory header — type aliases + Make() for the
// dret::pdl::DocListIdxPDL family.
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

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"
#include "dret/pdl/set_codecs.h"
#include "dret/pdl/storage_policy.h"

#include "../axes.h"
#include "../enum_traits.h"
#include "dgcda.h"
#include "gcda.h"
#include "slp_ns.h"

namespace bench::factories::pdl {

inline constexpr std::uint8_t kWidth = dret::Alphabet<>::int_width;

// GetDoc backing aliases. SLP / SLP-NS / DSLP are parameterised on the grammar
// type so PDL can pick any GCDA-family backend (the same TSLP defaults the RMQ
// and GCDA/DGCDA/SLP-NS factories use, so the typed cache files are shared).
template <typename TStorage>
using GetDocsDA   = dret::pdl::PDLGetDocsDA<TStorage, kWidth>;

template <typename TStorage, typename TSLP = gcda::SLP_Light>
using GetDocsSLP  = dret::pdl::PDLGetDocsSLP<TStorage, kWidth, TSLP>;

template <typename TStorage, typename TSLP = slp_ns::BareSLP_IV>
using GetDocsSLP_NS = dret::pdl::PDLGetDocsSLP_NS<TStorage, kWidth, TSLP>;

template <typename TStorage, typename TDSLP = dgcda::SLP_Default>
using GetDocsDSLP = dret::pdl::PDLGetDocsDSLP<TStorage, kWidth, TDSLP>;

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

template <typename T>
struct Tag { using type = T; };

// Build one PDL index for a fixed codec template, dispatching the get-doc
// backend (and its grammar-variant sub-axis) to the right TGetDocs at compile
// time. SLP fans across the GCDA TSLP variants, SLP-NS across the bare-SLP
// containers, DSLP across the DGCDA variants; DA ignores all three.
template <template <typename, typename> class TCodec, typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
MakeForCodec(TStorage t_storage, dret::Config& t_config,
             uint32_t t_block_size, float t_storing_factor,
             dret::pdl::StoragePolicy lib_policy,
             bench::axes::GetDocEnum t_get_doc,
             bench::axes::GCDASLPVariant t_gcda_slp,
             bench::axes::BareSLPVariant t_bare_slp,
             bench::axes::DGCDASLPVariant t_dgcda_slp) {
  using namespace bench::axes;
  std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t> result;

  auto build = [&](auto type_tag) {
    using TGetDocs = typename decltype(type_tag)::type;
    using TIndex = TCodec<TStorage, TGetDocs>;
    auto idx = std::make_shared<TIndex>(t_storage, t_block_size, t_storing_factor, lib_policy);
    construct(*idx, t_config);
    idx->load(t_config);
    result = {idx, sdsl::size_in_bytes(*idx)};
  };

  switch (t_get_doc) {
    case GetDocEnum::SLP:
      switch (t_gcda_slp) {
        case GCDASLPVariant::CompactBP:    build(Tag<GetDocsSLP<TStorage, gcda::SLP_CompactBP>>{}); break;
        case GCDASLPVariant::CompactLOUDS: build(Tag<GetDocsSLP<TStorage, gcda::SLP_CompactLOUDS>>{}); break;
        case GCDASLPVariant::Combined:     build(Tag<GetDocsSLP<TStorage, gcda::SLP_Combined>>{}); break;
        case GCDASLPVariant::Light:
        default:                           build(Tag<GetDocsSLP<TStorage, gcda::SLP_Light>>{}); break;
      }
      break;
    case GetDocEnum::SLP_NS:
      switch (t_bare_slp) {
        case BareSLPVariant::Raw: build(Tag<GetDocsSLP_NS<TStorage, slp_ns::BareSLP_Raw>>{}); break;
        case BareSLPVariant::DV:  build(Tag<GetDocsSLP_NS<TStorage, slp_ns::BareSLP_DV>>{});  break;
        case BareSLPVariant::VV:  build(Tag<GetDocsSLP_NS<TStorage, slp_ns::BareSLP_VV>>{});  break;
        case BareSLPVariant::IV:
        default:                  build(Tag<GetDocsSLP_NS<TStorage, slp_ns::BareSLP_IV>>{});  break;
      }
      break;
    case GetDocEnum::DSLP:
      switch (t_dgcda_slp) {
        case DGCDASLPVariant::OTF: build(Tag<GetDocsDSLP<TStorage, dgcda::SLP_OTF>>{}); break;
        case DGCDASLPVariant::CRL: build(Tag<GetDocsDSLP<TStorage, dgcda::SLP_CRL>>{}); break;
        case DGCDASLPVariant::EV:  build(Tag<GetDocsDSLP<TStorage, dgcda::SLP_EV>>{});  break;
        case DGCDASLPVariant::DV:  build(Tag<GetDocsDSLP<TStorage, dgcda::SLP_DV>>{});  break;
        case DGCDASLPVariant::VV:  build(Tag<GetDocsDSLP<TStorage, dgcda::SLP_VV>>{});  break;
        case DGCDASLPVariant::Default:
        default:                   build(Tag<GetDocsDSLP<TStorage, dgcda::SLP_Default>>{}); break;
      }
      break;
    case GetDocEnum::DA:
    default:
      build(Tag<GetDocsDA<TStorage>>{});
      break;
  }
  return result;
}

}  // namespace detail

// Factory-path entry. Dispatches codec, then the get-doc backend + its grammar
// variant; storage policy is threaded into the constructor.
template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage,
     dret::Config& t_config,
     uint32_t t_block_size,
     float t_storing_factor,
     bench::axes::PDLVariant t_codec,
     bench::axes::GetDocEnum t_get_doc,
     bench::axes::PDLStoragePolicy t_policy,
     bench::axes::GCDASLPVariant t_gcda_slp = bench::axes::GCDASLPVariant::Light,
     bench::axes::BareSLPVariant t_bare_slp = bench::axes::BareSLPVariant::IV,
     bench::axes::DGCDASLPVariant t_dgcda_slp = bench::axes::DGCDASLPVariant::Default) {
  using namespace bench::axes;
  const auto lib_policy = toPDLStoragePolicy(t_policy);
  switch (t_codec) {
    case PDLVariant::RP:
      return detail::MakeForCodec<IdxRP>(t_storage, t_config, t_block_size, t_storing_factor,
                                         lib_policy, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp);
    case PDLVariant::BC:
      return detail::MakeForCodec<IdxBC>(t_storage, t_config, t_block_size, t_storing_factor,
                                         lib_policy, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp);
    case PDLVariant::Plain:
    default:
      return detail::MakeForCodec<IdxPlain>(t_storage, t_config, t_block_size, t_storing_factor,
                                            lib_policy, t_get_doc, t_gcda_slp, t_bare_slp, t_dgcda_slp);
  }
}

}  // namespace bench::factories::pdl
