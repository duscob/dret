//
// GCDA factory header: the grammar and document-set types of the
// dret::gcda::DocListIdxGCDA family. Its cells are in gcda.cpp (and, over a
// document-array backend, gcda_backend.cpp); the RMQ and PDL factories reuse
// the grammar types, so their cache files are shared.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include <grammar/sampled_slp.h>
#include <grammar/slp.h>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_gcda_backend.h"
#include "dret/backend/da.h"
#include "dret/backend/rlcsa.h"
#include "dret/backend/sa_phi.h"
#include "dret/slp/combined_slp_with_unit_cover.h"
#include "dret/slp/compact_bp_slp.h"
#include "dret/slp/compact_louds_slp.h"

#include "../axes.h"

namespace bench::factories::gcda {

// TSLP type choices for the GCDA family. Storage-independent.
using SLP_Light     = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                          grammar::SampledSLP<>,
                                          grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
using SLP_CompactBP   = grammar::CompactBPSLP<>;
using SLP_CompactLOUDS = grammar::CompactLOUDSSLP<>;
// Bit-compressed (int_vector) rules, lengths, and leaves — matches light/compact.
// The all-default grammar::CombinedSLPWithUnitCover<> would store these as
// uncompressed std::vector<uint32_t> (~8x larger grammar base).
using SLP_Combined        = grammar::CombinedSLPWithUnitCover<
    grammar::CombinedSLP<grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
                         grammar::SampledSLP<>,
                         sdsl::int_vector<>>>;


// Document-set containers (GCDASetsCodec): Re-Pair compressed, the family
// default, or plain -- the sorted lists, bit-packed.
using Sets_RP = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                  true,
                                  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>;
using Sets_Plain = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>;

// Storage-parameterised typed-index template alias. The default TSLP matches
// the dret::gcda::DocListIdxGCDA default (SLP_Light = grammar::LightSLP<...>).
template <typename TStorage, typename TSLP = SLP_Light, typename TSets = Sets_RP>
using Idx = dret::gcda::DocListIdxGCDA<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP,
    TSets>;

// GCDA with the counts of its lists (computes document frequencies).
template <typename TStorage, typename TSLP = SLP_Light, typename TSets = Sets_RP>
using IdxCounts = dret::gcda::DocFreqIdxGCDA<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP,
    TSets>;

}  // namespace bench::factories::gcda
