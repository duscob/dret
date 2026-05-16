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

namespace bench::factories::pdl {

// GetDoc backing aliases — DA-backed, GCDA-LightSLP-backed, DGCDA-backed.
template <typename TStorage>
using GetDocsDA   = dret::pdl::PDLGetDocsDA<TStorage, dret::Alphabet<>::int_width>;

template <typename TStorage>
using GetDocsSLP  = dret::pdl::PDLGetDocsSLP<TStorage>;

template <typename TStorage>
using GetDocsDSLP = dret::pdl::PDLGetDocsDSLP<TStorage>;

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

// Factory-path entry. The codec × get-doc cross is dispatched via the
// type-tag-lambda pattern; storage policy is threaded into the constructor.
template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage,
     dret::Config& t_config,
     uint32_t t_block_size,
     float t_storing_factor,
     bench::axes::PDLVariant t_codec,
     bench::axes::GetDocEnum t_get_doc,
     bench::axes::PDLStoragePolicy t_policy) {
  using namespace bench::axes;
  const auto lib_policy = toPDLStoragePolicy(t_policy);
  std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t> result;

  auto build = [&](auto type_tag) {
    using TIndex = typename decltype(type_tag)::type;
    auto idx = std::make_shared<TIndex>(t_storage, t_block_size, t_storing_factor, lib_policy);
    construct(*idx, t_config);
    idx->load(t_config);
    result = {idx, sdsl::size_in_bytes(*idx)};
  };

  struct Plain_DA   { using type = IdxPlain<TStorage>; };
  struct Plain_SLP  { using type = IdxPlain<TStorage, GetDocsSLP<TStorage>>; };
  struct Plain_DSLP { using type = IdxPlain<TStorage, GetDocsDSLP<TStorage>>; };
  struct RP_DA      { using type = IdxRP<TStorage>; };
  struct RP_SLP     { using type = IdxRP<TStorage, GetDocsSLP<TStorage>>; };
  struct RP_DSLP    { using type = IdxRP<TStorage, GetDocsDSLP<TStorage>>; };
  struct BC_DA      { using type = IdxBC<TStorage>; };
  struct BC_SLP     { using type = IdxBC<TStorage, GetDocsSLP<TStorage>>; };
  struct BC_DSLP    { using type = IdxBC<TStorage, GetDocsDSLP<TStorage>>; };

  switch (t_codec) {
    case PDLVariant::Plain:
      switch (t_get_doc) {
        case GetDocEnum::SLP:  build(Plain_SLP{});  break;
        case GetDocEnum::DSLP: build(Plain_DSLP{}); break;
        case GetDocEnum::DA:
        default:               build(Plain_DA{});   break;
      }
      break;
    case PDLVariant::RP:
      switch (t_get_doc) {
        case GetDocEnum::SLP:  build(RP_SLP{});  break;
        case GetDocEnum::DSLP: build(RP_DSLP{}); break;
        case GetDocEnum::DA:
        default:               build(RP_DA{});   break;
      }
      break;
    case PDLVariant::BC:
      switch (t_get_doc) {
        case GetDocEnum::SLP:  build(BC_SLP{});  break;
        case GetDocEnum::DSLP: build(BC_DSLP{}); break;
        case GetDocEnum::DA:
        default:               build(BC_DA{});   break;
      }
      break;
  }
  return result;
}

}  // namespace bench::factories::pdl
