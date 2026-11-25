//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The PDL cells: codec x get-doc backend x its grammar variant x policy x block size x storing
// factor. The indexes are made in pdl_<codec>.cpp (pdl_codecs.h); this file only enumerates and
// names.
//

#include <string>
#include <vector>

#include "dret/config.h"
#include "dret/pdl/tree_core.h"

#include "../enum_traits.h"
#include "cells.h"
#include "names.h"
#include "pdl_codecs.h"

namespace bench::factories {

namespace {

using namespace bench::axes;
using namespace bench::names;

MakeIndex MakeFor(PDLVariant t_codec, const pdl::Params& t_params) {
  switch (t_codec) {
    case PDLVariant::RP:
      return pdl::MakeRP(t_params);
    case PDLVariant::BC:
      return pdl::MakeBC(t_params);
    case PDLVariant::Plain:
    default:
      return pdl::MakePlain(t_params);
  }
}

// The selection and the lists of the cell, in every codec; the tree (of a block size) is shared.
std::string RebuildPrefix(std::uint32_t t_bs, float t_sf, PDLStoragePolicy t_policy) {
  return dret::pdl::PdlSelectionPrefix(dret::kDefaultKeys.keys, t_bs, t_sf, toPDLStoragePolicy(t_policy)) + "pdl_";
}

}  // namespace

std::vector<Cell> Cells(const spec::PDLSweep& t_sweep) {
  std::vector<Cell> cells;
  for (auto codec : t_sweep.codec) {
    for (auto gd : t_sweep.get_doc) {
      // Every (policy, block size, storing factor) of one backend variant, its keys in t_variant_kv.
      auto emit = [&](const KV& t_variant_kv,
                      GCDASLPVariant t_gslp,
                      BareSLPVariant t_bslp,
                      DGCDASLPVariant t_dslp,
                      std::uint32_t t_backend_bs = 1024) {
        for (auto policy : t_sweep.policy)
          for (auto bs : t_sweep.block_size)
            for (auto sf : t_sweep.storing_factor) {
              KV kv{{"codec", CodecValue(codec)}, {"get-doc", GetDocValue(gd)}};
              for (const auto& p : t_variant_kv) kv.push_back(p);
              kv.push_back({"policy", PolicyValue(policy)});
              kv.push_back({"block-size", IntStr(bs)});
              kv.push_back({"storing-factor", SfStr(sf)});
              const pdl::Params params{bs, sf, policy, gd, t_gslp, t_bslp, t_dslp, t_backend_bs};
              Cell cell{.name = KeyedName("PDL", kv), .make = MakeFor(codec, params)};
              cell.rebuild_prefixes = {RebuildPrefix(bs, sf, policy)};
              cell.block_size = bs;
              cell.storing_factor = sf;
              cells.push_back(std::move(cell));
            }
      };
      const auto light = GCDASLPVariant::Light;
      const auto iv = BareSLPVariant::IV;
      const auto full_iv = DGCDASLPVariant::Default;
      switch (gd) {
        case GetDocEnum::SLP:
          // The backend block size joins the name only when it is not the default 1024.
          for (auto v : t_sweep.gcda_slp)
            for (auto bbs : t_sweep.backend_block_size) {
              KV vkv{{"slp", SlpPlainValue(v)}};
              if (bbs != 1024) vkv.push_back({"backend-block-size", IntStr(bbs)});
              emit(vkv, v, iv, full_iv, bbs);
            }
          break;
        case GetDocEnum::SLP_NS:
          // A differential backend samples every backend-block-size positions, always in its name: a
          // name without it means the spacing was PDL's block size, as before the backend had its own.
          for (auto v : t_sweep.bare_slp) {
            if (IsBareDiff(v)) {
              for (auto bbs : t_sweep.backend_block_size) {
                auto vkv = BareKeys(v);
                vkv.push_back({"backend-block-size", IntStr(bbs)});
                emit(vkv, light, v, full_iv, bbs);
              }
            } else {
              emit(BareKeys(v), light, v, full_iv);
            }
          }
          break;
        case GetDocEnum::DSLP:
          for (auto v : t_sweep.dgcda_slp) {
            auto [span, cont] = DiffParts(v);
            emit({{"span-length", span}, {"slp-container", cont}}, light, iv, v);
          }
          break;
        case GetDocEnum::SAPhiSR:
          // Named by the sr-index subsample rate; built as SA-Phi-R (see factories/pdl.h).
          for (auto s : t_sweep.sa_sampling) emit({{"sa-sampling", IntStr(s)}}, light, iv, full_iv);
          break;
        case GetDocEnum::SAPhiR:
        case GetDocEnum::RLCSA:
        case GetDocEnum::DA:
        default:
          emit({}, light, iv, full_iv);
          break;
      }
    }
  }
  return cells;
}

}  // namespace bench::factories
