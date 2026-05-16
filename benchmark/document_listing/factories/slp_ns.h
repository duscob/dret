//
// SLP-NS factory header — type aliases + Make() for the non-sampled
// dret::DocListIdxSLP family (Phase C). Also exposes the BareSLP_* type
// aliases that the RMQ-NS path in factories/rmq.h reuses to share cache
// files via grammar::SLP<TVars, TLens> type-hashing.
//

#pragma once

#include <cstddef>
#include <memory>
#include <utility>

#include <sdsl/dac_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/vlc_vector.hpp>

#include <grammar/slp.h>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_slp.h"

#include "../axes.h"

namespace bench::factories::slp_ns {

// Bare-SLP container choices. enc_vector<> is excluded: rule pairs and span
// lengths are non-monotonic. Default matches dret::DocListIdxSLP<>::TSLP.
using BareSLP_Default = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
using BareSLP_Raw     = grammar::SLP<>;
using BareSLP_DV      = grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>;
using BareSLP_VV      = grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>;

// Storage-parameterised typed-index template alias. Default TSLP matches
// dret::DocListIdxSLP<>::TSLP (BareSLP_Default).
template <typename TStorage, typename TSLP = BareSLP_Default>
using Idx = dret::DocListIdxSLP<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP>;

template <typename TStorage>
std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t>
Make(TStorage t_storage,
     dret::Config& t_config,
     bench::axes::BareSLPVariant t_slp) {
  std::pair<std::shared_ptr<dret::DocListIndex<>>, std::size_t> result;

  auto build = [&](auto type_tag) {
    using TIndex = typename decltype(type_tag)::type;
    auto idx = std::make_shared<TIndex>(t_storage);
    construct(*idx, t_config);
    idx->load(t_config);
    result = {idx, sdsl::size_in_bytes(*idx)};
  };

  struct T_Default { using type = Idx<TStorage, BareSLP_Default>; };
  struct T_Raw     { using type = Idx<TStorage, BareSLP_Raw>; };
  struct T_DV      { using type = Idx<TStorage, BareSLP_DV>; };
  struct T_VV      { using type = Idx<TStorage, BareSLP_VV>; };

  switch (t_slp) {
    case bench::axes::BareSLPVariant::Raw: build(T_Raw{}); break;
    case bench::axes::BareSLPVariant::DV:  build(T_DV{});  break;
    case bench::axes::BareSLPVariant::VV:  build(T_VV{});  break;
    case bench::axes::BareSLPVariant::Default:
    default:                               build(T_Default{}); break;
  }
  return result;
}

}  // namespace bench::factories::slp_ns
