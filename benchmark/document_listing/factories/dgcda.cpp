//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The cells of GCDA over a differential grammar: grammar variant x document-set codec x block size x
// storing factor.
//

#include <string>
#include <vector>

#include "cells.h"
#include "dgcda.h"
#include "names.h"

namespace bench::factories {

namespace {

using namespace bench::axes;

template <typename TSets>
MakeIndex Make(DGCDASLPVariant t_slp, std::uint32_t t_bs, float t_sf) {
  switch (t_slp) {
    case DGCDASLPVariant::OTF:
      return MakeOf<dgcda::Idx<Storage, dgcda::SLP_OTF, TSets>>(t_bs, t_sf);
    case DGCDASLPVariant::CRL:
      return MakeOf<dgcda::Idx<Storage, dgcda::SLP_CRL, TSets>>(t_bs, t_sf);
    case DGCDASLPVariant::EV:
      return MakeOf<dgcda::Idx<Storage, dgcda::SLP_EV, TSets>>(t_bs, t_sf);
    case DGCDASLPVariant::DV:
      return MakeOf<dgcda::Idx<Storage, dgcda::SLP_DV, TSets>>(t_bs, t_sf);
    case DGCDASLPVariant::VV:
      return MakeOf<dgcda::Idx<Storage, dgcda::SLP_VV, TSets>>(t_bs, t_sf);
    case DGCDASLPVariant::Default:
    default:
      return MakeOf<dgcda::Idx<Storage, dgcda::SLP_Default, TSets>>(t_bs, t_sf);
  }
}

}  // namespace

std::vector<Cell> Cells(const spec::DGCDASweep& t_sweep) {
  using namespace bench::names;
  const auto& keys = dret::kDefaultKeys.keys;
  std::vector<Cell> cells;
  for (auto sets : t_sweep.sets)
    for (auto slp : t_sweep.tslp)
      for (auto bs : t_sweep.block_size)
        for (auto sf : t_sweep.storing_factor) {
          auto [span, cont] = DiffParts(slp);
          Cell cell{.name = KeyedName("GCDA",
                                      {{"slp", "diff"},
                                       {"span-length", span},
                                       {"slp-container", cont},
                                       {"tree", "sampled"},
                                       {"sets", SetsValue(sets)},
                                       {"block-size", IntStr(bs)},
                                       {"storing-factor", SfStr(sf)}}),
                    .make = sets == GCDASetsCodec::Plain ? Make<dgcda::Sets_Plain>(slp, bs, sf)
                                                         : Make<dgcda::Sets_RP>(slp, bs, sf)};
          // The cell's components and the samples of its differential grammar, taken every bs positions.
          cell.rebuild_prefixes = {
              dret::PrefixedKey(keys, dret::conf::kBlkSf, "da_", bs, sf),
              dret::PrefixedKey(keys, dret::conf::kSpc, dret::KeyName(keys, dret::conf::kDaDiffSamples), bs)};
          cell.block_size = bs;
          cell.storing_factor = sf;
          cells.push_back(std::move(cell));
        }
  return cells;
}

}  // namespace bench::factories
