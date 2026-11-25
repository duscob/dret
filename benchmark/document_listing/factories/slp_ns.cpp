//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The Brute-GCDA cells: brute force over a grammar without a tree (tree=none). A bare grammar (plain,
// or differential, sampled every block-size positions), or a sampled GCDA grammar without its
// document lists, at the (block size, storing factor) of the GCDA cell it reads.
//

#include <string>
#include <vector>

#include "cells.h"
#include "names.h"
#include "slp_ns.h"

namespace bench::factories {

namespace {

using namespace bench::axes;

// A bare grammar ignores the storing factor; a plain one the spacing too.
MakeIndex MakeBare(BareSLPVariant t_slp, std::uint32_t t_spacing) {
  constexpr float kSf = dret::kNolistsStoringFactor;
  switch (t_slp) {
    case BareSLPVariant::Raw:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_Raw>>(t_spacing, kSf);
    case BareSLPVariant::DV:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_DV>>(t_spacing, kSf);
    case BareSLPVariant::VV:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_VV>>(t_spacing, kSf);
    case BareSLPVariant::Diff:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_Diff>>(t_spacing, kSf);
    case BareSLPVariant::DiffEV:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_DiffEV>>(t_spacing, kSf);
    case BareSLPVariant::DiffDV:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_DiffDV>>(t_spacing, kSf);
    case BareSLPVariant::DiffVV:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_DiffVV>>(t_spacing, kSf);
    case BareSLPVariant::IV:
    default:
      return MakeOf<slp_ns::Idx<Storage, slp_ns::BareSLP_IV>>(t_spacing, kSf);
  }
}

MakeIndex MakeSampled(GCDASLPVariant t_slp, std::uint32_t t_bs, float t_sf) {
  switch (t_slp) {
    case GCDASLPVariant::Combined:
      return MakeOf<slp_ns::Idx<Storage, gcda::SLP_Combined>>(t_bs, t_sf);
    case GCDASLPVariant::CompactBP:
      return MakeOf<slp_ns::Idx<Storage, gcda::SLP_CompactBP>>(t_bs, t_sf);
    case GCDASLPVariant::CompactLOUDS:
      return MakeOf<slp_ns::Idx<Storage, gcda::SLP_CompactLOUDS>>(t_bs, t_sf);
    case GCDASLPVariant::Light:
    default:
      return MakeOf<slp_ns::Idx<Storage, gcda::SLP_Light>>(t_bs, t_sf);
  }
}

}  // namespace

std::vector<Cell> Cells(const spec::SLPNSSweep& t_sweep) {
  using namespace bench::names;
  const auto& keys = dret::kDefaultKeys.keys;
  std::vector<Cell> cells;
  for (auto slp : t_sweep.tslp) {
    // A differential grammar sweeps its sample spacing; a plain one has none and is one cell. The
    // spacing joins the name only when it is not 512.
    const bool diff = IsBareDiff(slp);
    const auto spacings = diff ? t_sweep.block_size : std::vector<std::uint32_t>{512};
    for (auto bs : spacings) {
      const auto name = diff
                            ? KeyedName("GCDA",
                                        {{"slp", "diff"},
                                         {"slp-container", DiffContainerValue(slp)},
                                         {"tree", "none"},
                                         {"block-size", bs == 512 ? std::string() : IntStr(bs)}})
                            : KeyedName("GCDA", {{"slp", "bare"}, {"slp-container", BareValue(slp)}, {"tree", "none"}});
      Cell cell{.name = name, .make = MakeBare(slp, bs)};
      // The plain grammar is the CNF grammar itself; a differential one owns its samples.
      cell.rebuild_prefixes = {
          dret::KeyName(keys, dret::conf::kDaCnfGrammar) + "_",
          dret::PrefixedKey(keys, dret::conf::kSpc, dret::KeyName(keys, dret::conf::kDaDiffSamples), bs)};
      if (diff) cell.block_size = bs;
      cells.push_back(std::move(cell));
    }
  }
  // Sampled GCDA grammars without their lists. --rebuild leaves them alone: they are GCDA's own
  // components.
  for (auto slp : t_sweep.gcda_slp)
    for (auto bs : t_sweep.block_size)
      for (auto sf : t_sweep.storing_factor) {
        Cell cell{.name = KeyedName("GCDA",
                                    {{"slp", SlpPlainValue(slp)},
                                     {"tree", "none"},
                                     {"block-size", IntStr(bs)},
                                     {"storing-factor", SfStr(sf)}}),
                  .make = MakeSampled(slp, bs, sf)};
        cell.block_size = bs;
        cell.storing_factor = sf;
        cells.push_back(std::move(cell));
      }
  return cells;
}

}  // namespace bench::factories
