//
// Brute-baseline factory header — type aliases for the dret::DocListIdxBrute
// family wrapping an r-index or sr-index.
//
// Brute baselines bypass the construct-then-load auto-build pattern: their
// underlying sri::RIndex / sri::SrIndexValidArea caches are managed
// externally (in the query factory: via the lazy rIndex() / srIndex()
// accessors; in bm_build_items: via DocListIdxBrute::construct() which
// builds the locate-index from the raw data file). The factory.h facade
// special-cases brute and does not call construct() on the index itself —
// only load(). bm_build_items uses DocListIdxBrute<> directly with
// BM_ConstructBruteIdx<> which exercises construct(idx, config).
//

#pragma once

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "dret/doc_list/doc_list_brute.h"

namespace bench::factories::brute {

// r-index-backed brute baseline.
template <typename TStorage>
using IdxR = dret::DocListIdxBrute<
    TStorage,
    dret::Alphabet<>,
    sri::RIndex<TStorage>>;

// sr-index-backed brute baseline. The sampling parameter is threaded
// through the index constructor as a runtime argument.
template <typename TStorage>
using IdxSr = dret::DocListIdxBrute<
    TStorage,
    dret::Alphabet<>,
    sri::SrIndexValidArea<TStorage>>;

// The same shape, named for construct mode. DocListIdxBrute's own default
// locate index is sri::SrIdxGeneric<SrIndexValidArea<...>, 16>, which binds
// the subsample rate at compile time; construct mode has to cover whatever
// sampling_size values the spec lists, so it names SrIndexValidArea directly
// and passes the rate to the constructor. That is all SrIdxGeneric does
// internally -- it exists only to expose a one-argument (storage) ctor.
template <typename TStorage>
using IdxSrConstruct = IdxSr<TStorage>;

}  // namespace bench::factories::brute
