//
// bm_doc_list — JSON-spec-driven doc-list benchmark binary.
//
// Reads a JSON sweep spec (see spec.h for the schema), expands each family
// block's Cartesian product, and registers one Google Benchmark per cell.
// The 'mode' field selects between query-time and construction-time
// benchmarking; the two modes share the spec, the dataset / workload
// settings, and the per-family sweep blocks.
//
// Usage:
//   bm_doc_list --spec=sweep.json [--benchmark_out_format=csv --benchmark_out=path.csv]
//

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <benchmark/benchmark.h>
#include <gflags/gflags.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "../bm_size_counters.h"

#include "../tool/definitions.h"

#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/size_report.h"

#include "axes.h"
#include "enum_traits.h"
#include "factories/brute.h"
#include "factories/dgcda.h"
#include "factories/gcda.h"
#include "factories/pdl.h"
#include "factories/rmq.h"
#include "factories/slp_ns.h"
#include "factory.h"
#include "spec.h"

DEFINE_string(spec, "", "Path to the JSON sweep spec (MANDATORY).");

namespace {

namespace fac = bench::factories;
using bench::axes::BareSLPVariant;
using bench::axes::DGCDASLPVariant;
using bench::axes::GCDASLPVariant;
using bench::axes::GetDocEnum;
using bench::axes::PDLStoragePolicy;
using bench::axes::PDLVariant;
using bench::axes::EnumTraits;
using GS = dret::GenericStorage;
using ExternalGS = std::reference_wrapper<sri::GenericStorage>;

//~~~~~~~  Query-mode result containers (vector of doc ids) ~~~~~~~

struct DocListResult {
  std::vector<std::size_t> docs;
  void operator()(std::size_t d) { docs.push_back(d); }
  void operator()() {
    std::sort(docs.begin(), docs.end());
    docs.erase(std::unique(docs.begin(), docs.end()), docs.end());
  }
};

//~~~~~~~  Counter helpers (inlined to avoid pulling in bm_query.h's BM_WarmUp) ~~~~~~~

void SetupQueryCounters(benchmark::State& t_state) {
  t_state.counters["Size(bytes)"] = 0;
  t_state.counters["Bits_x_Symbol"] = 0;
  t_state.counters["Patterns"] = 0;
  t_state.counters["Time_x_Pattern"] = 0;
}

//~~~~~~~  Query benchmark kernel ~~~~~~~

auto BM_Query = [](benchmark::State& t_state,
                   Factory<>* t_factory,
                   Factory<>::Config t_config,
                   const std::vector<std::string>* t_patterns,
                   std::size_t t_seq_size) {
  auto [idx, idx_size] = t_factory->Make(t_config);

  for (auto _ : t_state) {
    for (const auto& pattern : *t_patterns) {
      DocListResult result;
      idx->Search(pattern, std::ref(result));
      result();
      benchmark::DoNotOptimize(result);
    }
  }

  SetupQueryCounters(t_state);
  t_state.counters["Size(bytes)"] = idx_size;
  t_state.counters["Bits_x_Symbol"] = idx_size * 8.0 / t_seq_size;
  t_state.counters["Patterns"] = t_patterns->size();
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_patterns->size(),
      benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);

  if (idx) {
    const auto sizes = idx->GetSizeReport();
    appendCounters(t_state, sizes);
    t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  }
};

//~~~~~~~  Construct benchmark kernels ~~~~~~~
//
// Construct kernels build a fresh, owning index and time the construct() call.
// The benchmark instance owns a GenericStorage so the typed index isn't tied
// to the factory's external-storage cache.

void SetupConstructCounters(benchmark::State& t_state, const dret::Config& t_config,
                             std::uint32_t bs, float sf) {
  using namespace sri::conf;
  sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
  t_state.counters["n"] = buf.size();
  if (bs > 0) t_state.counters["bs"] = bs;
  if (sf > 0) t_state.counters["sf"] = sf;
}

template <typename TIndex>
void BM_ConstructBrute(benchmark::State& t_state, dret::Config t_config) {
  TIndex index;
  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }
  SetupConstructCounters(t_state, t_config, 0, 0);
}

template <typename TIndex>
void BM_ConstructGCDAFamily(benchmark::State& t_state, dret::Config t_config) {
  const auto bs = static_cast<std::uint32_t>(t_state.range(0));
  const auto sf = static_cast<float>(t_state.range(1));
  GS storage;
  TIndex index(storage, bs, sf);
  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }
  SetupConstructCounters(t_state, t_config, bs, sf);
  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
}

template <typename TIndex>
void BM_ConstructSLPNS(benchmark::State& t_state, dret::Config t_config) {
  GS storage;
  TIndex index(storage);
  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }
  SetupConstructCounters(t_state, t_config, 0, 0);
  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
}

template <typename TIndex, typename TCore>
void BM_ConstructRMQ(benchmark::State& t_state, dret::Config t_config) {
  const auto bs = static_cast<std::uint32_t>(t_state.range(0));
  const auto sf = static_cast<float>(t_state.range(1));
  GS storage;
  TCore core(storage, bs, sf);
  TIndex index(storage, core);
  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }
  SetupConstructCounters(t_state, t_config, bs, sf);
  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
}

template <typename TIndex, typename TCore>
void BM_ConstructRMQ_NS(benchmark::State& t_state, dret::Config t_config) {
  GS storage;
  TCore core(storage);
  TIndex index(storage, core);
  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }
  SetupConstructCounters(t_state, t_config, 0, 0);
  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
}

template <typename TIndex>
void BM_ConstructPDL(benchmark::State& t_state, dret::Config t_config,
                     dret::pdl::StoragePolicy t_policy) {
  const auto bs = static_cast<std::uint32_t>(t_state.range(0));
  const auto sf = static_cast<float>(t_state.range(1));
  GS storage;
  TIndex index(storage, bs, sf, t_policy);
  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    construct(index, t_config);
    sdsl::memory_monitor::stop();
  }
  SetupConstructCounters(t_state, t_config, bs, sf);
  index.load(t_config);
  const auto sizes = index.GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
}

//~~~~~~~  Helpers ~~~~~~~

std::string NameSuffix(GCDASLPVariant v) {
  return (v == GCDASLPVariant::Light)
      ? std::string{}
      : std::string("-") + EnumTraits<GCDASLPVariant>::Name(v);
}

std::string NameSuffix(DGCDASLPVariant v) {
  return (v == DGCDASLPVariant::Default)
      ? std::string{}
      : std::string("-") + EnumTraits<DGCDASLPVariant>::Name(v);
}

std::string NameSuffix(BareSLPVariant v) {
  return (v == BareSLPVariant::IV)
      ? std::string{}
      : std::string("-") + EnumTraits<BareSLPVariant>::Name(v);
}

const char* CoreName(fac::rmq::CoreKind c) {
  switch (c) {
    case fac::rmq::CoreKind::SADA:  return "SADA";
    case fac::rmq::CoreKind::ILCP:  return "ILCP";
    case fac::rmq::CoreKind::CILCP: return "CILCP";
  }
  return "UNKNOWN";
}

std::string BsSfSuffix(std::uint32_t bs, float sf) {
  return "-bs" + std::to_string(bs) + "-sf" + std::to_string(static_cast<int>(sf));
}

std::vector<std::int64_t> ArgsRange(const std::vector<std::uint32_t>& bs,
                                     const std::vector<float>& sf) {
  // Google Benchmark's ArgsProduct takes vector<vector<int64_t>>.
  return {};
}

//~~~~~~~  Query-mode register functions (one per family) ~~~~~~~

void RegisterQueryGCDA(const bench::spec::GCDASweep& sw, Factory<>& factory,
                       const std::vector<std::string>* patterns, std::size_t seq_size) {
  for (auto tslp : sw.tslp) {
    for (auto bs : sw.block_size) {
      for (auto sf : sw.storing_factor) {
        const auto name = "DocListGCDA" + NameSuffix(tslp) + BsSfSuffix(bs, sf);
        Factory<>::Config cfg{};
        cfg.index_t = Factory<>::IndexEnum::GCDA;
        cfg.block_size = bs;
        cfg.storing_factor = sf;
        cfg.gcda_slp = tslp;
        benchmark::RegisterBenchmark(name, BM_Query, &factory, cfg, patterns, seq_size);
      }
    }
  }
}

void RegisterQueryDGCDA(const bench::spec::DGCDASweep& sw, Factory<>& factory,
                        const std::vector<std::string>* patterns, std::size_t seq_size) {
  for (auto tslp : sw.tslp) {
    for (auto bs : sw.block_size) {
      for (auto sf : sw.storing_factor) {
        const auto name = "DocListDGCDA" + NameSuffix(tslp) + BsSfSuffix(bs, sf);
        Factory<>::Config cfg{};
        cfg.index_t = Factory<>::IndexEnum::DGCDA;
        cfg.block_size = bs;
        cfg.storing_factor = sf;
        cfg.dgcda_slp = tslp;
        benchmark::RegisterBenchmark(name, BM_Query, &factory, cfg, patterns, seq_size);
      }
    }
  }
}

void RegisterQuerySLPNS(const bench::spec::SLPNSSweep& sw, Factory<>& factory,
                        const std::vector<std::string>* patterns, std::size_t seq_size) {
  for (auto tslp : sw.tslp) {
    const auto name = "DocListSLP-NS" + NameSuffix(tslp);
    Factory<>::Config cfg{};
    cfg.index_t = Factory<>::IndexEnum::SLP_NS;
    cfg.bare_slp = tslp;
    benchmark::RegisterBenchmark(name, BM_Query, &factory, cfg, patterns, seq_size);
  }
}

void RegisterQueryRMQ(const bench::spec::RMQSweep& sw, Factory<>& factory,
                      const std::vector<std::string>* patterns, std::size_t seq_size) {
  for (auto core : sw.core) {
    Factory<>::IndexEnum core_idx{};
    switch (core) {
      case fac::rmq::CoreKind::SADA:  core_idx = Factory<>::IndexEnum::SADA;  break;
      case fac::rmq::CoreKind::ILCP:  core_idx = Factory<>::IndexEnum::ILCP;  break;
      case fac::rmq::CoreKind::CILCP: core_idx = Factory<>::IndexEnum::CILCP; break;
    }
    for (auto gd : sw.get_doc) {
      const bool needs_bs_sf = (gd == GetDocEnum::SLP) || (gd == GetDocEnum::DSLP);
      auto bs_list = needs_bs_sf ? sw.block_size : std::vector<std::uint32_t>{0};
      auto sf_list = needs_bs_sf ? sw.storing_factor : std::vector<float>{0.0f};
      // Inner SLP axis: gcda_slp for get_doc=slp; bare_slp for slp_ns; ignored otherwise.
      auto inner_label = [&](auto variant) -> std::string {
        if (gd == GetDocEnum::SLP)    return NameSuffix(std::get<0>(variant));
        if (gd == GetDocEnum::SLP_NS) return NameSuffix(std::get<1>(variant));
        return {};
      };
      auto inner_vec = [&]() -> std::vector<std::pair<GCDASLPVariant, BareSLPVariant>> {
        if (gd == GetDocEnum::SLP) {
          std::vector<std::pair<GCDASLPVariant, BareSLPVariant>> v;
          for (auto x : sw.gcda_slp) v.push_back({x, BareSLPVariant::IV});
          return v;
        }
        if (gd == GetDocEnum::SLP_NS) {
          std::vector<std::pair<GCDASLPVariant, BareSLPVariant>> v;
          for (auto x : sw.bare_slp) v.push_back({GCDASLPVariant::Light, x});
          return v;
        }
        return {{GCDASLPVariant::Light, BareSLPVariant::IV}};
      }();
      for (auto bs : bs_list) {
        for (auto sf : sf_list) {
          for (auto pair : inner_vec) {
            std::string name = std::string(CoreName(core)) + "-"
                             + EnumTraits<GetDocEnum>::Name(gd)
                             + inner_label(pair);
            if (needs_bs_sf) name += BsSfSuffix(bs, sf);
            Factory<>::Config cfg{};
            cfg.index_t = core_idx;
            cfg.get_doc = gd;
            cfg.gcda_slp = pair.first;
            cfg.bare_slp = pair.second;
            cfg.block_size = needs_bs_sf ? bs : 512;
            cfg.storing_factor = needs_bs_sf ? sf : 4.0f;
            benchmark::RegisterBenchmark(name, BM_Query, &factory, cfg, patterns, seq_size);
          }
        }
      }
    }
  }
}

void RegisterQueryPDL(const bench::spec::PDLSweep& sw, Factory<>& factory,
                      const std::vector<std::string>* patterns, std::size_t seq_size) {
  for (auto codec : sw.codec) {
    for (auto gd : sw.get_doc) {
      for (auto policy : sw.policy) {
        for (auto bs : sw.block_size) {
          for (auto sf : sw.storing_factor) {
            std::string name = std::string("DocListPDL-")
                             + EnumTraits<PDLVariant>::Name(codec) + "-"
                             + EnumTraits<GetDocEnum>::Name(gd) + "-"
                             + EnumTraits<PDLStoragePolicy>::Name(policy)
                             + BsSfSuffix(bs, sf);
            Factory<>::Config cfg{};
            cfg.index_t = Factory<>::IndexEnum::PDL;
            cfg.block_size = bs;
            cfg.storing_factor = sf;
            cfg.get_doc = gd;
            cfg.pdl_variant = codec;
            cfg.pdl_storage_policy = policy;
            benchmark::RegisterBenchmark(name, BM_Query, &factory, cfg, patterns, seq_size);
          }
        }
      }
    }
  }
}

void RegisterQueryBrute(const bench::spec::BruteSweep& sw, Factory<>& factory,
                        const std::vector<std::string>* patterns, std::size_t seq_size) {
  using bench::spec::BruteSweep;
  for (auto kind : sw.kind) {
    if (kind == BruteSweep::Kind::RIndex) {
      Factory<>::Config cfg{};
      cfg.index_t = Factory<>::IndexEnum::BRUTE_R_INDEX;
      benchmark::RegisterBenchmark("DocListBrute-RIndex", BM_Query,
                                    &factory, cfg, patterns, seq_size);
    } else {
      for (auto s : sw.sampling_size) {
        Factory<>::Config cfg{};
        cfg.index_t = Factory<>::IndexEnum::BRUTE_SR_INDEX;
        cfg.sampling_size = s;
        const auto name = "DocListBrute-SRIndex-s" + std::to_string(s);
        benchmark::RegisterBenchmark(name, BM_Query, &factory, cfg, patterns, seq_size);
      }
    }
  }
}

//~~~~~~~  Construct-mode register functions (one per family) ~~~~~~~

template <typename TIndex>
void RegisterOneConstructGCDA(const std::string& name, dret::Config& config,
                               const std::vector<std::uint32_t>& bs,
                               const std::vector<float>& sf) {
  std::vector<std::int64_t> bs_i64(bs.begin(), bs.end());
  std::vector<std::int64_t> sf_i64(sf.begin(), sf.end());
  benchmark::RegisterBenchmark(name, BM_ConstructGCDAFamily<TIndex>, config)
      ->ArgsProduct({bs_i64, sf_i64});
}

void RegisterConstructGCDA(const bench::spec::GCDASweep& sw, dret::Config& config) {
  using namespace fac::gcda;
  for (auto tslp : sw.tslp) {
    const auto name = "DocListGCDA" + NameSuffix(tslp);
    switch (tslp) {
      case GCDASLPVariant::CompactBP:
        RegisterOneConstructGCDA<Idx<GS, SLP_CompactBP>>(name, config, sw.block_size, sw.storing_factor); break;
      case GCDASLPVariant::CompactLOUDS:
        RegisterOneConstructGCDA<Idx<GS, SLP_CompactLOUDS>>(name, config, sw.block_size, sw.storing_factor); break;
      case GCDASLPVariant::Combined:
        RegisterOneConstructGCDA<Idx<GS, SLP_Combined>>(name, config, sw.block_size, sw.storing_factor); break;
      case GCDASLPVariant::Light:
      default:
        RegisterOneConstructGCDA<Idx<GS>>(name, config, sw.block_size, sw.storing_factor); break;
    }
  }
}

void RegisterConstructDGCDA(const bench::spec::DGCDASweep& sw, dret::Config& config) {
  using namespace fac::dgcda;
  for (auto tslp : sw.tslp) {
    const auto name = "DocListDGCDA" + NameSuffix(tslp);
    switch (tslp) {
      case DGCDASLPVariant::OTF: RegisterOneConstructGCDA<Idx<GS, SLP_OTF>>(name, config, sw.block_size, sw.storing_factor); break;
      case DGCDASLPVariant::CRL: RegisterOneConstructGCDA<Idx<GS, SLP_CRL>>(name, config, sw.block_size, sw.storing_factor); break;
      case DGCDASLPVariant::EV:  RegisterOneConstructGCDA<Idx<GS, SLP_EV>>(name, config, sw.block_size, sw.storing_factor); break;
      case DGCDASLPVariant::DV:  RegisterOneConstructGCDA<Idx<GS, SLP_DV>>(name, config, sw.block_size, sw.storing_factor); break;
      case DGCDASLPVariant::VV:  RegisterOneConstructGCDA<Idx<GS, SLP_VV>>(name, config, sw.block_size, sw.storing_factor); break;
      case DGCDASLPVariant::Default:
      default:                   RegisterOneConstructGCDA<Idx<GS>>(name, config, sw.block_size, sw.storing_factor); break;
    }
  }
}

void RegisterConstructSLPNS(const bench::spec::SLPNSSweep& sw, dret::Config& config) {
  using namespace fac::slp_ns;
  for (auto tslp : sw.tslp) {
    const auto name = "DocListSLP-NS" + NameSuffix(tslp);
    switch (tslp) {
      case BareSLPVariant::Raw: benchmark::RegisterBenchmark(name, BM_ConstructSLPNS<Idx<GS, BareSLP_Raw>>, config); break;
      case BareSLPVariant::DV:  benchmark::RegisterBenchmark(name, BM_ConstructSLPNS<Idx<GS, BareSLP_DV>>, config); break;
      case BareSLPVariant::VV:  benchmark::RegisterBenchmark(name, BM_ConstructSLPNS<Idx<GS, BareSLP_VV>>, config); break;
      case BareSLPVariant::IV:
      default:                  benchmark::RegisterBenchmark(name, BM_ConstructSLPNS<Idx<GS>>, config); break;
    }
  }
}

// RMQ construct dispatch: per (core, get_doc, [gcda_slp | bare_slp]) tuple.
// Helper sub-templates pick the typed TIndex / TCore for the cell.
template <template <typename, typename> class TCoreT, typename TGetDoc, bool NeedsBsSf>
void RegisterRMQOneTriple(const std::string& core_name, const std::string& suffix,
                          dret::Config& config,
                          const std::vector<std::uint32_t>& bs,
                          const std::vector<float>& sf) {
  using TCore = TCoreT<GS, TGetDoc>;
  using TIndex = fac::rmq::Idx<GS, TCore>;
  const auto name = core_name + suffix;
  if constexpr (NeedsBsSf) {
    std::vector<std::int64_t> bs_i64(bs.begin(), bs.end());
    std::vector<std::int64_t> sf_i64(sf.begin(), sf.end());
    benchmark::RegisterBenchmark(name, BM_ConstructRMQ<TIndex, TCore>, config)
        ->ArgsProduct({bs_i64, sf_i64});
  } else {
    benchmark::RegisterBenchmark(name, BM_ConstructRMQ_NS<TIndex, TCore>, config);
  }
}

template <template <typename, typename> class TCoreT>
void RegisterRMQOneCore(const std::string& core_name,
                         const bench::spec::RMQSweep& sw, dret::Config& config) {
  using namespace fac::rmq;
  for (auto gd : sw.get_doc) {
    const bool needs_bs_sf = (gd == GetDocEnum::SLP) || (gd == GetDocEnum::DSLP);
    switch (gd) {
      case GetDocEnum::DA: {
        // DA-backed: one entry per core, no bs/sf.
        using TCore = TCoreT<GS, GetDocDA<GS>>;
        using TIndex = Idx<GS, TCore>;
        benchmark::RegisterBenchmark(core_name + "-DA",
            BM_ConstructBrute<TIndex>, config);
        break;
      }
      case GetDocEnum::SLP: {
        for (auto tslp : sw.gcda_slp) {
          const auto suffix = std::string("-SLP") + NameSuffix(tslp);
          switch (tslp) {
            case GCDASLPVariant::CompactBP:    RegisterRMQOneTriple<TCoreT, GetDocSLP<GS, SLP_CompactBP>, true>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
            case GCDASLPVariant::CompactLOUDS: RegisterRMQOneTriple<TCoreT, GetDocSLP<GS, SLP_CompactLOUDS>, true>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
            case GCDASLPVariant::Combined:     RegisterRMQOneTriple<TCoreT, GetDocSLP<GS, SLP_Combined>, true>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
            case GCDASLPVariant::Light:
            default:                           RegisterRMQOneTriple<TCoreT, GetDocSLP<GS>, true>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
          }
        }
        break;
      }
      case GetDocEnum::SLP_NS: {
        for (auto tslp : sw.bare_slp) {
          const auto suffix = std::string("-SLP-NS") + NameSuffix(tslp);
          switch (tslp) {
            case BareSLPVariant::Raw: RegisterRMQOneTriple<TCoreT, GetDocSLP_NS<GS, BareSLP_Raw>, false>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
            case BareSLPVariant::DV:  RegisterRMQOneTriple<TCoreT, GetDocSLP_NS<GS, BareSLP_DV>, false>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
            case BareSLPVariant::VV:  RegisterRMQOneTriple<TCoreT, GetDocSLP_NS<GS, BareSLP_VV>, false>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
            case BareSLPVariant::IV:
            default:                  RegisterRMQOneTriple<TCoreT, GetDocSLP_NS<GS>, false>(core_name, suffix, config, sw.block_size, sw.storing_factor); break;
          }
        }
        break;
      }
      case GetDocEnum::DSLP: {
        RegisterRMQOneTriple<TCoreT, GetDocDSLP<GS>, true>(core_name, "-DSLP", config, sw.block_size, sw.storing_factor);
        break;
      }
    }
  }
}

void RegisterConstructRMQ(const bench::spec::RMQSweep& sw, dret::Config& config) {
  using namespace fac::rmq;
  for (auto core : sw.core) {
    switch (core) {
      case CoreKind::SADA:  RegisterRMQOneCore<SadaCore>("DocListSADA", sw, config);  break;
      case CoreKind::ILCP:  RegisterRMQOneCore<IlcpCore>("DocListILCP", sw, config);  break;
      case CoreKind::CILCP: RegisterRMQOneCore<CilcpCore>("DocListCILCP", sw, config); break;
    }
  }
}

void RegisterConstructPDL(const bench::spec::PDLSweep& sw, dret::Config& config) {
  using namespace fac::pdl;
  std::vector<std::int64_t> bs_i64(sw.block_size.begin(), sw.block_size.end());
  std::vector<std::int64_t> sf_i64(sw.storing_factor.begin(), sw.storing_factor.end());
  for (auto codec : sw.codec) {
    for (auto gd : sw.get_doc) {
      for (auto policy : sw.policy) {
        const auto lib_policy = bench::axes::toPDLStoragePolicy(policy);
        const auto name = std::string("DocListPDL-")
                        + EnumTraits<PDLVariant>::Name(codec) + "-"
                        + EnumTraits<GetDocEnum>::Name(gd) + "-"
                        + EnumTraits<PDLStoragePolicy>::Name(policy);
        auto reg = [&]<typename TIndex>() {
          benchmark::RegisterBenchmark(name, BM_ConstructPDL<TIndex>, config, lib_policy)
              ->ArgsProduct({bs_i64, sf_i64});
        };
        switch (codec) {
          case PDLVariant::Plain:
            switch (gd) {
              case GetDocEnum::DA:   reg.template operator()<IdxPlain<GS>>(); break;
              case GetDocEnum::SLP:  reg.template operator()<IdxPlain<GS, GetDocsSLP<GS>>>(); break;
              case GetDocEnum::DSLP: reg.template operator()<IdxPlain<GS, GetDocsDSLP<GS>>>(); break;
              case GetDocEnum::SLP_NS: break;
            }
            break;
          case PDLVariant::RP:
            switch (gd) {
              case GetDocEnum::DA:   reg.template operator()<IdxRP<GS>>(); break;
              case GetDocEnum::SLP:  reg.template operator()<IdxRP<GS, GetDocsSLP<GS>>>(); break;
              case GetDocEnum::DSLP: reg.template operator()<IdxRP<GS, GetDocsDSLP<GS>>>(); break;
              case GetDocEnum::SLP_NS: break;
            }
            break;
          case PDLVariant::BC:
            switch (gd) {
              case GetDocEnum::DA:   reg.template operator()<IdxBC<GS>>(); break;
              case GetDocEnum::SLP:  reg.template operator()<IdxBC<GS, GetDocsSLP<GS>>>(); break;
              case GetDocEnum::DSLP: reg.template operator()<IdxBC<GS, GetDocsDSLP<GS>>>(); break;
              case GetDocEnum::SLP_NS: break;
            }
            break;
        }
      }
    }
  }
}

void RegisterConstructBrute(const bench::spec::BruteSweep& sw, dret::Config& config,
                            const std::string& data_path) {
  using bench::spec::BruteSweep;
  for (auto kind : sw.kind) {
    if (kind == BruteSweep::Kind::RIndex) {
      benchmark::RegisterBenchmark("DocListIdxBrute",
          BM_ConstructBrute<dret::DocListIdxBrute<>>, config);
    }
    // sr-index construct path is not exposed today; skip.
  }
}

}  // namespace

int main(int argc, char** argv) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_spec.empty()) {
    std::cerr << "ERROR: --spec=<path.json> is required" << std::endl;
    return 1;
  }

  bench::spec::Spec spec;
  try {
    spec = bench::spec::ParseSpecFile(FLAGS_spec);
  } catch (const std::exception& e) {
    std::cerr << "ERROR parsing spec: " << e.what() << std::endl;
    return 2;
  }

  // Build the dret::Config from the dataset block.
  dret::Config config;
  std::string data_path;
  if (spec.mode == bench::spec::Mode::Query) {
    config = dret::Config(spec.dataset.name, spec.dataset.dir,
                          sri::SDSL_LIBDIVSUFSORT, true);
    data_path = spec.dataset.dir + "/" + spec.dataset.name;
  } else {
    data_path = spec.dataset.dir + "/" + spec.dataset.name;
    config = dret::Config(data_path,
                          std::filesystem::current_path(),
                          sri::toSAAlgo(spec.dataset.sa_algo),
                          false,
                          spec.dataset.data_width,
                          spec.dataset.doc_delim);
  }

  // Read patterns up-front for query mode (kept alive through main).
  std::vector<std::string> patterns;
  if (spec.mode == bench::spec::Mode::Query) {
    std::ifstream pf(spec.dataset.patterns, std::ios::binary);
    if (!pf) {
      std::cerr << "ERROR: cannot open patterns file: " << spec.dataset.patterns << std::endl;
      return 3;
    }
    std::string line;
    while (std::getline(pf, line)) {
      if (!line.empty()) patterns.push_back(line);
    }
  }

  std::size_t seq_size = 0;
  std::unique_ptr<Factory<>> factory;
  if (spec.mode == bench::spec::Mode::Query) {
    factory = std::make_unique<Factory<>>(config);
    seq_size = factory->SequenceSize();
  }

  // Dispatch each sweep block per mode.
  for (const auto& fs : spec.sweep) {
    if (spec.mode == bench::spec::Mode::Query) {
      std::visit([&]<typename T>(const T& sw) {
        if constexpr (std::is_same_v<T, bench::spec::GCDASweep>) RegisterQueryGCDA(sw, *factory, &patterns, seq_size);
        else if constexpr (std::is_same_v<T, bench::spec::DGCDASweep>) RegisterQueryDGCDA(sw, *factory, &patterns, seq_size);
        else if constexpr (std::is_same_v<T, bench::spec::SLPNSSweep>) RegisterQuerySLPNS(sw, *factory, &patterns, seq_size);
        else if constexpr (std::is_same_v<T, bench::spec::RMQSweep>) RegisterQueryRMQ(sw, *factory, &patterns, seq_size);
        else if constexpr (std::is_same_v<T, bench::spec::PDLSweep>) RegisterQueryPDL(sw, *factory, &patterns, seq_size);
        else if constexpr (std::is_same_v<T, bench::spec::BruteSweep>) RegisterQueryBrute(sw, *factory, &patterns, seq_size);
      }, fs);
    } else {
      std::visit([&]<typename T>(const T& sw) {
        if constexpr (std::is_same_v<T, bench::spec::GCDASweep>) RegisterConstructGCDA(sw, config);
        else if constexpr (std::is_same_v<T, bench::spec::DGCDASweep>) RegisterConstructDGCDA(sw, config);
        else if constexpr (std::is_same_v<T, bench::spec::SLPNSSweep>) RegisterConstructSLPNS(sw, config);
        else if constexpr (std::is_same_v<T, bench::spec::RMQSweep>) RegisterConstructRMQ(sw, config);
        else if constexpr (std::is_same_v<T, bench::spec::PDLSweep>) RegisterConstructPDL(sw, config);
        else if constexpr (std::is_same_v<T, bench::spec::BruteSweep>) RegisterConstructBrute(sw, config, data_path);
      }, fs);
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  return 0;
}
