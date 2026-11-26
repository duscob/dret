//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The RMQ frequency cells (dret/doc_freq/rmq_core.h): a published core x a frequency oracle. With the
// rank oracle the backend is the wavelet tree over the DA; with the doc-ISA oracle it is the
// differential suffix array, and the block size is the sample spacing of both grammars. They compute
// document frequencies; a listing query gets the published core's answer over the same backend.
//

#include <stdexcept>
#include <string>
#include <vector>

#include "dret/backend/diff_sa.h"
#include "dret/backend/wavelet_da.h"
#include "dret/doc_freq/oracle.h"
#include "dret/doc_freq/rmq_core.h"

#include "cells.h"
#include "names.h"

namespace bench::factories {

namespace {

using namespace bench::axes;

constexpr auto kW = dret::Alphabet<>::int_width;
using Count = sri::RIndexCount<Storage, dret::Alphabet<>>;
using Cores = dret::freq::Cores<Storage, kW>;
using WT = dret::backend::WaveletDA<Storage, kW>;
using DSA = dret::backend::DiffSA<Storage, kW>;

template <template <typename> class TCoreOf>
MakeIndex Make(FreqOracle t_oracle, std::uint32_t t_bs) {
  if (t_oracle == FreqOracle::Rank)
    return MakeOf<dret::freq::DocFreqIdxRMQ<Storage, dret::Alphabet<>, Count, WT, TCoreOf, dret::freq::WtRankOracle<Storage, kW>>>();
  return MakeOf<dret::freq::DocFreqIdxRMQ<Storage, dret::Alphabet<>, Count, DSA, TCoreOf, dret::freq::DocIsaOracle<Storage, kW>>>(
      t_bs);
}

}  // namespace

std::vector<Cell> Cells(const spec::RMQFreqSweep& t_sweep) {
  using namespace bench::names;
  std::vector<Cell> cells;
  for (auto core : t_sweep.core) {
    if (core != CoreKind::SADA && core != CoreKind::ILCP && core != CoreKind::CILCP)
      throw std::invalid_argument("rmq_freq: core must be sada, ilcp or cilcp (the published cores)");
    for (auto oracle : t_sweep.oracle) {
      const bool rank = oracle == FreqOracle::Rank;
      // The block size names only the doc-ISA cells, the ones it changes.
      const auto bs_list = rank ? std::vector<std::uint32_t>{0} : t_sweep.block_size;
      for (auto bs : bs_list) {
        KV kv{{"get-doc", rank ? "wt-da" : "dsa"}, {"freq", rank ? "rank" : "doc-isa"}};
        if (!rank) kv.push_back({"block-size", IntStr(bs)});
        Cell cell{.name = KeyedName(CoreName(core), kv),
                  .make = core == CoreKind::SADA   ? Make<Cores::Sada>(oracle, bs)
                          : core == CoreKind::ILCP ? Make<Cores::Ilcp>(oracle, bs)
                                                   : Make<Cores::Cilcp>(oracle, bs)};
        // No rebuild: the reversed arrays, the wavelet tree and the grammars are shared by the cores.
        cells.push_back(std::move(cell));
      }
    }
  }
  return cells;
}

}  // namespace bench::factories
