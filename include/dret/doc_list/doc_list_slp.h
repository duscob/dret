//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// Brute-GCDA: brute force over a GCDA grammar without its document lists (brute_core.h over
// backend::BareGrammar). The grammar is the CNF grammar alone (a bare grammar::SLP: the range is split
// into a span cover, each variable expanded), the differential grammar (range decompression from the
// nearest sample), or a sampled GCDA grammar (cached, on-demand, compact) whose lists are not loaded
// (expansion through its sampled tree, as the RMQ/PDL GCDA backend does).
//
// Against GCDA it has no stored lists: the whole [sp, ep) is expanded, O((ep - sp) * height) plus the
// cover, and each document reported once by the listing bitmap (distinct_docs.h).
//
// The parameters are the sample spacing of a differential grammar, or the block size and storing factor
// of a sampled one; a plain grammar ignores both.
//

#pragma once

#include <grammar/slp.h>
#include <sdsl/int_vector.hpp>

#include "sr-index/r_index.h"

#include "dret/backend/grammar.h"
#include "dret/config.h"
#include "dret/doc_list/brute_core.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/slp/bare_grammar.h"

namespace dret {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>
using DocListIdxSLP =
    DocListIdx<TStorage,
               TAlphabet,
               TCountIdx,
               BruteCore<TStorage, TAlphabet::int_width, backend::BareGrammar<TStorage, TAlphabet::int_width, TSLP>>>;

}  // namespace dret
