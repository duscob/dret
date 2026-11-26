//
// DGCDA factory header: the differential grammar types of the
// dret::dgcda::DocListIdxDGCDA family (DocListIdxGCDA with TSLP defaulting to
// DifferentialLightSLP<>). Its cells are in dgcda.cpp.
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
// All variants store the SLP base rules/lengths bit-compressed (int_vector):
// SLP_Default via DifferentialLightSLP<>'s default TSLP; the OTF/CRL adapters
// over BasicSLP<int_vector>; EV/DV/VV via an explicit int_vector base.
using DSLP_Base = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
using SLP_Default = dret::DifferentialLightSLP<>;
using SLP_OTF     = dret::DifferentialLightSLP<
    dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<sdsl::int_vector<>>>>;
using SLP_CRL     = dret::DifferentialLightSLP<
    dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<sdsl::int_vector<>>>>;
using SLP_EV = dret::DifferentialLightSLP<DSLP_Base,
                                          grammar::SampledSLP<>,
                                          sdsl::enc_vector<>,
                                          sdsl::enc_vector<>,
                                          sdsl::enc_vector<>>;
using SLP_DV = dret::DifferentialLightSLP<DSLP_Base,
                                          grammar::SampledSLP<>,
                                          sdsl::dac_vector<>,
                                          sdsl::dac_vector<>,
                                          sdsl::dac_vector<>>;
using SLP_VV = dret::DifferentialLightSLP<DSLP_Base,
                                          grammar::SampledSLP<>,
                                          sdsl::vlc_vector<>,
                                          sdsl::vlc_vector<>,
                                          sdsl::vlc_vector<>>;

// Storage-parameterised typed-index template alias. Default TSLP matches
// dret::dgcda::DocListIdxDGCDA<>::TSLP (SLP_Default).

// Document-set containers (GCDASetsCodec): Re-Pair compressed, the family
// default, or plain -- the sorted lists, bit-packed.
using Sets_RP = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                  true,
                                  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
using Sets_Plain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;

template <typename TStorage, typename TSLP = SLP_Default, typename TSets = Sets_RP>
using Idx = dret::dgcda::DocListIdxDGCDA<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP,
    TSets>;

}  // namespace bench::factories::dgcda
