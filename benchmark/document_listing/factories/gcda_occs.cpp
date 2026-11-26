//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The GCDA cells with occurrence bitvectors (dret/doc_freq/gcda_occs_core.h), one per frequency
// oracle: the wavelet tree over the DA with its rank, or the differential suffix array with the
// per-document ISAs, whose grammars' sample spacing is the block size. They compute document
// frequencies; a listing query gets the documents.
//

#include <string>
#include <vector>

#include "sr-index/r_index.h"

#include "dret/backend/diff_sa.h"
#include "dret/backend/wavelet_da.h"
#include "dret/doc_freq/gcda_occs_core.h"
#include "dret/doc_freq/oracle.h"

#include "cells.h"
#include "names.h"

namespace bench::factories {

std::vector<Cell> Cells(const spec::GCDAOccsSweep& t_sweep) {
  using namespace bench::names;
  using namespace bench::axes;
  constexpr auto kW = dret::Alphabet<>::int_width;
  using Count = sri::RIndexCount<Storage, dret::Alphabet<>>;
  std::vector<Cell> cells;
  for (auto oracle : t_sweep.oracle) {
    const bool rank = oracle == FreqOracle::Rank;
    // The block size names only the doc-ISA cells, the ones it changes.
    const auto bs_list = rank ? std::vector<std::uint32_t>{0} : t_sweep.block_size;
    for (auto bs : bs_list) {
      KV kv{{"tree", "none"}, {"occs", "plain"}, {"get-doc", rank ? "wt-da" : "dsa"}, {"freq", rank ? "rank" : "doc-isa"}};
      if (!rank) kv.push_back({"block-size", IntStr(bs)});
      Cell cell{.name = KeyedName("GCDA", kv),
                .make = rank ? MakeOf<dret::freq::DocFreqIdxGCDAOccs<Storage, dret::Alphabet<>, Count,
                                                                      dret::backend::WaveletDA<Storage, kW>,
                                                                      dret::freq::WtRankOracle<Storage, kW>>>()
                             : MakeOf<dret::freq::DocFreqIdxGCDAOccs<Storage, dret::Alphabet<>, Count,
                                                                      dret::backend::DiffSA<Storage, kW>,
                                                                      dret::freq::DocIsaOracle<Storage, kW>>>(bs)};
      // No rebuild: the grammar, the bitvectors and the oracle structures are shared.
      cells.push_back(std::move(cell));
    }
  }
  return cells;
}

}  // namespace bench::factories
