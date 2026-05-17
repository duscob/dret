//
// DGCDA factory header — type aliases + Make() for the
// dret::dgcda::DocListIdxDGCDA family (the alias of DocListIdxGCDA with
// TSLP defaulting to DifferentialLightSLP<>).
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/vlc_vector.hpp>

#include <grammar/sampled_slp.h>
#include <grammar/slp.h>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/slp/basic_slp_span_length.h"
#include "dret/slp/differential_light_slp.h"

#include "../axes.h"

namespace bench::factories::dgcda {

// TSLP type choices for the DGCDA family. Six variants:
// - Default: standard DifferentialLightSLP<>
// - OTF: on-the-fly span lengths (no stored lengths)
// - CRL: cached root span lengths only
// - EV / DV / VV: vary the inner int-vector container of roots / span_sums /
//   samples (TSampleRootsPos keeps its default sdsl::enc_vector<>).
using SLP_Default = dret::DifferentialLightSLP<>;
using SLP_OTF     = dret::DifferentialLightSLP<
    dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>;
using SLP_CRL     = dret::DifferentialLightSLP<
    dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>;
using SLP_EV = dret::DifferentialLightSLP<grammar::SLP<>,
                                          grammar::SampledSLP<>,
                                          sdsl::enc_vector<>,
                                          sdsl::enc_vector<>,
                                          sdsl::enc_vector<>>;
using SLP_DV = dret::DifferentialLightSLP<grammar::SLP<>,
                                          grammar::SampledSLP<>,
                                          sdsl::dac_vector<>,
                                          sdsl::dac_vector<>,
                                          sdsl::dac_vector<>>;
using SLP_VV = dret::DifferentialLightSLP<grammar::SLP<>,
                                          grammar::SampledSLP<>,
                                          sdsl::vlc_vector<>,
                                          sdsl::vlc_vector<>,
                                          sdsl::vlc_vector<>>;

// Storage-parameterised typed-index template alias. Default TSLP matches
// dret::dgcda::DocListIdxDGCDA<>::TSLP (SLP_Default).
template <typename TStorage, typename TSLP = SLP_Default>
using Idx = dret::dgcda::DocListIdxDGCDA<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP>;

// Factory-path entry. Builds + loads a typed-index variant selected by the
// DGCDASLPVariant axis and returns the type-erased DocListIndex<> handle.
template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage,
     dret::Config& t_config,
     uint32_t t_block_size,
     float t_storing_factor,
     bench::axes::DGCDASLPVariant t_slp) {
  std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t> result;

  auto build = [&](auto type_tag) {
    using TIndex = typename decltype(type_tag)::type;
    auto idx = std::make_shared<TIndex>(t_storage, t_block_size, t_storing_factor);
    construct(*idx, t_config);
    idx->load(t_config);
    result = {idx, sdsl::size_in_bytes(*idx)};
  };

  struct T_Default { using type = Idx<TStorage, SLP_Default>; };
  struct T_OTF     { using type = Idx<TStorage, SLP_OTF>; };
  struct T_CRL     { using type = Idx<TStorage, SLP_CRL>; };
  struct T_EV      { using type = Idx<TStorage, SLP_EV>; };
  struct T_DV      { using type = Idx<TStorage, SLP_DV>; };
  struct T_VV      { using type = Idx<TStorage, SLP_VV>; };

  switch (t_slp) {
    case bench::axes::DGCDASLPVariant::OTF: build(T_OTF{}); break;
    case bench::axes::DGCDASLPVariant::CRL: build(T_CRL{}); break;
    case bench::axes::DGCDASLPVariant::EV:  build(T_EV{});  break;
    case bench::axes::DGCDASLPVariant::DV:  build(T_DV{});  break;
    case bench::axes::DGCDASLPVariant::VV:  build(T_VV{});  break;
    case bench::axes::DGCDASLPVariant::Default:
    default:                                build(T_Default{}); break;
  }
  return result;
}

}  // namespace bench::factories::dgcda
