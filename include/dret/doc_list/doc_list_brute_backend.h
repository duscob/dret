//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/19/26.
//
// Brute force over a backend (brute_core.h): find the suffix-array range [sp, ep) of the pattern and
// report the document of every position in it, each once. The backend is any document-array backend
// (backend/backend.h): the plain DA, SA-Phi, a GCDA grammar, the RLCSA.
//
// A backend that finds the range itself (the RLCSA) replaces the r-index, as in PDL: count_idx is then
// neither loaded, serialized nor counted in the size (doc_list_index.h).
//

#pragma once

#include "sr-index/r_index.h"

#include "dret/backend/da.h"
#include "dret/config.h"
#include "dret/doc_list/brute_core.h"
#include "dret/doc_list/doc_list_index.h"

namespace dret {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TGetDoc = backend::PlainDA<TStorage, TAlphabet::int_width>>
using DocListIdxBruteBackend =
    DocListIdx<TStorage, TAlphabet, TCountIdx, BruteCore<TStorage, TAlphabet::int_width, TGetDoc>>;

}  // namespace dret
