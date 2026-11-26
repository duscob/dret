//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The GCDA cells over its own grammar: grammar representation x document-set codec x block size x
// storing factor.
//

#include <string>
#include <vector>

#include "cells.h"
#include "gcda.h"
#include "names.h"

namespace bench::factories {

namespace {

using namespace bench::axes;

// TIdx: gcda::Idx, or gcda::IdxCounts for the cells with counts.
template <template <typename, typename, typename> class TIdx, typename TSets>
MakeIndex Make(GCDASLPVariant t_slp, std::uint32_t t_bs, float t_sf) {
  switch (t_slp) {
    case GCDASLPVariant::CompactBP:
      return MakeOf<TIdx<Storage, gcda::SLP_CompactBP, TSets>>(t_bs, t_sf);
    case GCDASLPVariant::CompactLOUDS:
      return MakeOf<TIdx<Storage, gcda::SLP_CompactLOUDS, TSets>>(t_bs, t_sf);
    case GCDASLPVariant::Combined:
      return MakeOf<TIdx<Storage, gcda::SLP_Combined, TSets>>(t_bs, t_sf);
    case GCDASLPVariant::Light:
    default:
      return MakeOf<TIdx<Storage, gcda::SLP_Light, TSets>>(t_bs, t_sf);
  }
}

MakeIndex Make(GCDASLPVariant t_slp, GCDASetsCodec t_sets, GCDACounts t_counts, std::uint32_t t_bs, float t_sf) {
  const bool plain = t_sets == GCDASetsCodec::Plain;
  if (t_counts == GCDACounts::Plain)
    return plain ? Make<gcda::IdxCounts, gcda::Sets_Plain>(t_slp, t_bs, t_sf)
                 : Make<gcda::IdxCounts, gcda::Sets_RP>(t_slp, t_bs, t_sf);
  return plain ? Make<gcda::Idx, gcda::Sets_Plain>(t_slp, t_bs, t_sf) : Make<gcda::Idx, gcda::Sets_RP>(t_slp, t_bs, t_sf);
}

}  // namespace

std::vector<Cell> Cells(const spec::GCDASweep& t_sweep) {
  using namespace bench::names;
  std::vector<Cell> cells;
  for (auto counts : t_sweep.counts)
    for (auto sets : t_sweep.sets)
      for (auto slp : t_sweep.tslp)
        for (auto bs : t_sweep.block_size)
          for (auto sf : t_sweep.storing_factor) {
            Cell cell{.name = KeyedName("GCDA",
                                        {{"slp", SlpPlainValue(slp)},
                                         {"tree", "sampled"},
                                         {"sets", SetsValue(sets)},
                                         {"counts", CountsValue(counts)},
                                         {"block-size", IntStr(bs)},
                                         {"storing-factor", SfStr(sf)}}),
                      .make = Make(slp, sets, counts, bs, sf)};
            // The cell's own components; the DA grammars and every shared structure stay.
            cell.rebuild_prefixes = {dret::PrefixedKey(dret::kDefaultKeys.keys, dret::conf::kBlkSf, "da_", bs, sf)};
            cell.block_size = bs;
            cell.storing_factor = sf;
            cells.push_back(std::move(cell));
          }
  return cells;
}

}  // namespace bench::factories
