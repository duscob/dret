//
// SLP-NS factory header: the bare grammar types of the dret::DocListIdxSLP
// family (Brute-GCDA). Its cells are in slp_ns.cpp; the RMQ and PDL factories
// reuse the BareSLP_* types, so their cache files are shared through
// grammar::SLP<TVars, TLens> type-hashing.
//

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <utility>

#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/vlc_vector.hpp>

#include <grammar/slp.h>

#include "sr-index/r_index.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_slp.h"

#include "../axes.h"
#include "gcda.h"

namespace bench::factories::slp_ns {

// Bare-SLP container choices. enc_vector<> is excluded: rule pairs and span
// lengths are non-monotonic. Default matches dret::DocListIdxSLP<>::TSLP.
using BareSLP_IV = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
using BareSLP_Raw     = grammar::SLP<>;
using BareSLP_DV      = grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>;
using BareSLP_VV      = grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>;
// bare-diff: base differential SLP, non-sampled (no GCChunks), sampled every
// spacing positions; its grammar is shared by every spacing. The container axis (roots/
// span_sums/samples) mirrors the DGCDA sampled-diff EV/DV/VV variants.
// Base SLP rules/lengths are bit-compressed (int_vector, = BareSLP_IV); the
// container suffix varies the per-field roots/span_sums/samples encoding.
using BareSLP_Diff    = dret::DifferentialSLP<BareSLP_IV>;  // int_vector (iv)
using BareSLP_DiffEV  = dret::DifferentialSLP<BareSLP_IV, sdsl::enc_vector<>, sdsl::enc_vector<>, sdsl::enc_vector<>>;
using BareSLP_DiffDV  = dret::DifferentialSLP<BareSLP_IV, sdsl::dac_vector<>, sdsl::dac_vector<>, sdsl::dac_vector<>>;
using BareSLP_DiffVV  = dret::DifferentialSLP<BareSLP_IV, sdsl::vlc_vector<>, sdsl::vlc_vector<>, sdsl::vlc_vector<>>;

// Storage-parameterised typed-index template alias. Default TSLP matches
// dret::DocListIdxSLP<>::TSLP (BareSLP_IV).
template <typename TStorage, typename TSLP = BareSLP_IV>
using Idx = dret::DocListIdxSLP<
    TStorage,
    dret::Alphabet<>,
    sri::RIndexCount<TStorage, dret::Alphabet<>>,
    TSLP>;

}  // namespace bench::factories::slp_ns
