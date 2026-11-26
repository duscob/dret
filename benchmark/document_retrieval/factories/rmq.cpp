//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The RMQ cells: core x get-doc backend x its grammar variant x the core's stored values. The
// indexes are made in rmq_<core>.cpp (rmq_cores.h); this file only enumerates and names.
//

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cells.h"
#include "names.h"
#include "rmq_cores.h"

namespace bench::factories {

namespace {

using namespace bench::axes;
using namespace bench::names;

MakeIndex MakeFor(CoreKind t_core, const rmq::Params& t_params) {
  switch (t_core) {
    case CoreKind::SADA:
      return rmq::MakeSada(t_params);
    case CoreKind::ILCP:
      return rmq::MakeIlcp(t_params);
    case CoreKind::CILCP:
      return rmq::MakeCilcp(t_params);
    case CoreKind::SADA_V:
      return rmq::MakeSadaValue(t_params);
    case CoreKind::ILCP_V:
      return rmq::MakeIlcpValue(t_params);
    case CoreKind::CILCP_V:
      return rmq::MakeCilcpValue(t_params);
  }
  throw std::invalid_argument("rmq: unknown core");
}

// The cache files a core's construct() writes, deleted by --rebuild. The grammar / DA backends are
// shared with the GCDA / DGCDA / brute cells and stay. rmq_n_doc is a one-element vector every RMQ
// build regenerates. The published (marking) cores own the shared structures; each -V core reuses
// them and owns only the array it adds (SADA-V reuses sada_rmq_; ILCP-V and CILCP-V reuse the rmq and
// run heads of ILCP and CILCP, which build the same partitions), so a --rebuild of a -V cell does not
// rebuild the partition: rebuild the published core for that.
std::vector<std::string> RebuildPrefixes(CoreKind t_core) {
  switch (t_core) {
    case CoreKind::SADA:
      return {"sada_rmq_", "rmq_n_doc_"};
    case CoreKind::ILCP:
      return {"ilcp_rmq_", "ilcp_run_heads_", "rmq_n_doc_"};
    case CoreKind::CILCP:
      return {"cilcp_rmq_", "cilcp_run_heads_", "rmq_n_doc_"};
    case CoreKind::SADA_V:
      return {"sada_prev_doc_", "rmq_n_doc_"};
    case CoreKind::ILCP_V:
      return {"ilcp_run_values_", "rmq_n_doc_"};
    case CoreKind::CILCP_V:
      return {"cilcp_run_values_", "rmq_n_doc_"};
  }
  return {};
}

}  // namespace

std::vector<Cell> Cells(const spec::RMQSweep& t_sweep) {
  std::vector<Cell> cells;
  for (auto core : t_sweep.core) {
    // ILCP-V / CILCP-V fan out across their stored run values, SADA-V across prev_doc; the published
    // cores have neither.
    const bool run_values_axis = core == CoreKind::ILCP_V || core == CoreKind::CILCP_V;
    const bool prev_doc_axis = core == CoreKind::SADA_V;
    const auto rv_list = run_values_axis ? t_sweep.run_values : std::vector<RunValuesVariant>{RunValuesVariant::DV};
    const auto pd_list = prev_doc_axis ? t_sweep.prev_doc : std::vector<PrevDocVariant>{PrevDocVariant::IV};

    for (auto gd : t_sweep.get_doc) {
      if (gd == GetDocEnum::RLCSA)
        throw std::invalid_argument("rmq: get_doc=rlcsa is not supported (see factories/rmq.h)");
      // A sampled grammar backend has a block size and a storing factor; a bare differential one a
      // block size (its sample spacing); the others neither.
      const bool sampled = gd == GetDocEnum::SLP || gd == GetDocEnum::DSLP;
      const bool any_bare_diff =
          gd == GetDocEnum::SLP_NS && std::any_of(t_sweep.bare_slp.begin(), t_sweep.bare_slp.end(), IsBareDiff);
      const auto bs_list = sampled || any_bare_diff ? t_sweep.block_size : std::vector<std::uint32_t>{0};
      const auto sf_list = sampled ? t_sweep.storing_factor : std::vector<float>{0.0f};
      // The grammar variant of the backend: gcda_slp for get_doc=slp, bare_slp for slp_ns.
      std::vector<std::pair<GCDASLPVariant, BareSLPVariant>> variants;
      if (gd == GetDocEnum::SLP)
        for (auto v : t_sweep.gcda_slp) variants.push_back({v, BareSLPVariant::IV});
      else if (gd == GetDocEnum::SLP_NS)
        for (auto v : t_sweep.bare_slp) variants.push_back({GCDASLPVariant::Light, v});
      else
        variants.push_back({GCDASLPVariant::Light, BareSLPVariant::IV});
      const auto dslp_list =
          gd == GetDocEnum::DSLP ? t_sweep.dgcda_slp : std::vector<DGCDASLPVariant>{DGCDASLPVariant::Default};
      // The sr-index subsample rate names SA-Phi-SR cells (built as SA-Phi-R, see factories/rmq.h).
      const auto sa_list = gd == GetDocEnum::SAPhiSR ? t_sweep.sa_sampling : std::vector<std::size_t>{0};

      for (auto bs : bs_list)
        for (auto sf : sf_list)
          for (auto [gslp, bslp] : variants) {
            const bool bare_diff = gd == GetDocEnum::SLP_NS && IsBareDiff(bslp);
            // A bare plain grammar ignores the block size: one cell.
            if (gd == GetDocEnum::SLP_NS && !bare_diff && bs != bs_list.front()) continue;
            for (auto dslp : dslp_list)
              for (auto rv : rv_list)
                for (auto pd : pd_list)
                  for (auto sa : sa_list) {
                    KV kv{{"get-doc", GetDocValue(gd)}};
                    if (gd == GetDocEnum::SLP) {
                      kv.push_back({"slp", SlpPlainValue(gslp)});
                    } else if (gd == GetDocEnum::SLP_NS) {
                      for (auto& p : BareKeys(bslp)) kv.push_back(p);
                    } else if (gd == GetDocEnum::DSLP) {
                      auto [span, cont] = DiffParts(dslp);
                      kv.push_back({"span-length", span});
                      kv.push_back({"slp-container", cont});
                    } else if (gd == GetDocEnum::SAPhiSR) {
                      kv.push_back({"sa-sampling", IntStr(sa)});
                    }
                    if (sampled || bare_diff) kv.push_back({"block-size", IntStr(bs)});
                    if (sampled) kv.push_back({"storing-factor", SfStr(sf)});
                    if (prev_doc_axis) kv.push_back({"prev-doc", PrevDocValue(pd)});
                    if (run_values_axis) kv.push_back({"run-values", RunValuesValue(rv)});

                    const std::uint32_t cell_bs = sampled || bare_diff ? bs : 512;
                    const float cell_sf = sampled ? sf : 4.0f;
                    const rmq::Params params{cell_bs, cell_sf, gd, gslp, bslp, dslp, rv, pd};
                    Cell cell{.name = KeyedName(CoreName(core), kv), .make = MakeFor(core, params)};
                    cell.rebuild_prefixes = RebuildPrefixes(core);
                    if (sampled || bare_diff) cell.block_size = bs;
                    if (sampled) cell.storing_factor = sf;
                    cells.push_back(std::move(cell));
                  }
          }
    }
  }
  return cells;
}

}  // namespace bench::factories
