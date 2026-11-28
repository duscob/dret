//
// bm_doc_list — JSON-spec-driven doc-list benchmark binary.
//
// Reads a JSON sweep spec (see spec.h for the schema), turns each family block
// into its cells (factories/cells.h: one per index, named FAMILY[k=v|...]), and
// registers one Google Benchmark per cell. The 'mode' field selects query-time
// or construction-time benchmarking; a cell has the same name in both modes.
//
// Usage:
//   bm_doc_list --spec=sweep.json [--benchmark_out_format=csv --benchmark_out=path.csv]
//

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#endif

#include <benchmark/benchmark.h>
#include <gflags/gflags.h>

#include <sdsl/config.hpp>
#include <sdsl/io.hpp>
#include <sdsl/memory_management.hpp>

#include "sr-index/config.h"

#include "../bm_size_counters.h"

#include "dret/config.h"
#include "dret/size_report.h"

#include "factories/cells.h"
#include "factory.h"
#include "spec.h"

#ifdef DRET_DOC_LIST_PROFILE
#include "dret/doc_list/search_profile.h"
#endif

DEFINE_string(spec, "", "Path to the JSON sweep spec (MANDATORY).");
DEFINE_string(print_results_dir, "",
              "If set, write each query index's per-pattern sorted doc-id lists to "
              "<dir>/<sanitized-name>.txt — all indexes must produce identical files "
              "(cross-index correctness verification; brute is the ground truth).");
DEFINE_string(per_pattern_dir, "",
              "If set, time every pattern of each query benchmark and write "
              "<dir>/<sanitized-name>.call<k>.tsv, one per call of the benchmark "
              "function (pattern index, length, occs, ndocs, minimum "
              "time in ns), plus Time_x_Pattern_p50/_p90/_p99 counters. Occs come from "
              "the r-index, so BruteRI's cache files must exist.");


namespace {

using bench::factories::Cell;

//~~~~~~~  Cache helpers for the construct-mode --rebuild plumbing ~~~~~~~

namespace cache_clean {

// Delete every regular file under cache_dir named <file_prefix>...<file_suffix>.
// SDSL writes cache files as "<key>_<basename>.sdsl", so the cell's key prefix
// anchors the front and "_<basename>.sdsl" anchors the back. The shared
// artefacts (Text / SA / DocEnds / DA / LCP) don't start with a per-cell key
// prefix, so they survive.
// Returns the number of files removed, so the caller can detect a prefix that
// matches nothing — a silently ineffective rebuild is indistinguishable from a
// fast one by timing alone.
std::size_t DeleteByPrefix(const std::filesystem::path& cache_dir,
                           const std::string& file_prefix,
                           const std::string& file_suffix) {
  namespace fs = std::filesystem;
  std::error_code ec;
  std::size_t removed = 0;
  if (!fs::is_directory(cache_dir, ec)) return removed;
  for (const auto& entry : fs::directory_iterator(cache_dir, ec)) {
    if (ec) break;
    const auto name = entry.path().filename().string();
    if (name.size() < file_prefix.size() + file_suffix.size()) continue;
    if (name.rfind(file_prefix, 0) != 0) continue;  // starts_with
    if (name.compare(name.size() - file_suffix.size(), file_suffix.size(),
                     file_suffix) != 0) continue;  // ends_with
    if (fs::remove(entry.path(), ec)) ++removed;
  }
  return removed;
}

// Build a rebuild closure for one construct-mode cell. Empty when disabled
// or when no prefixes are given. The closure deletes every cache file named
// <key_prefix>..._<basename>.sdsl for any prefix in the list.
// Multi-prefix support: RMQ cells need to wipe both the per-core RMQ cache
// (e.g. "sada_rmq_") and the shared rmq_n_doc helper (a one-element file
// rebuilt almost for free), but never the SLP cache (shared with GCDA).
std::function<void()> MakeHook(bool enable,
                               std::filesystem::path cache_dir,
                               std::string basename,
                               std::vector<std::string> key_prefixes) {
  if (!enable || key_prefixes.empty()) return {};
  const std::string file_suffix = "_" + basename + ".sdsl";
  return [cache_dir = std::move(cache_dir), file_suffix,
          key_prefixes = std::move(key_prefixes),
          first_call = std::make_shared<bool>(true)]() {
    std::size_t removed = 0;
    for (const auto& kp : key_prefixes) {
      removed += DeleteByPrefix(cache_dir, kp, file_suffix);
    }
    // A rebuild that deletes nothing yields a warm-cache measurement while
    // looking like a clean one. Report it on the first pass rather than let
    // the cell masquerade as a cold construct.
    if (*first_call && removed == 0) {
      std::cerr << "WARNING: rebuild deleted no cache files for prefixes [";
      for (std::size_t i = 0; i < key_prefixes.size(); ++i) {
        std::cerr << (i ? ", " : "") << key_prefixes[i];
      }
      std::cerr << "] with suffix '" << file_suffix << "' in " << cache_dir
                << " — timing will be cache-warm.\n";
    }
    *first_call = false;
  };
}

}  // namespace cache_clean


// Benchmark name -> file name: keep [A-Za-z0-9.-], map everything else to '_'.
std::string SanitizeName(std::string t_name) {
  for (char& c : t_name)
    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-' || c == '.')) c = '_';
  return t_name;
}

// Warmth check: stderr-warn once per cell when --rebuild is off AND the
// first construct() iteration returned under 1 ms (almost certainly a
// cache-warm no-op, not a real measurement).
void WarnIfWarm(const benchmark::State& t_state,
                std::int64_t first_iter_ns, bool rebuild_enabled) {
  constexpr std::int64_t kThresholdNs = 1'000'000;  // 1 ms
  if (rebuild_enabled || first_iter_ns < 0 || first_iter_ns >= kThresholdNs) return;
  std::cerr << "WARNING: " << t_state.name()
            << " construct() returned in " << first_iter_ns
            << " ns on the first iteration — cache was likely warm. "
            << "Set workload.rebuild=true in the spec for clean timing.\n";
}

// Memory-monitor trace + size sidecar emit, gated on workload.memory_trace.
// Files land in the working directory keyed by the sanitized benchmark name.
void WriteMemoryTrace(const benchmark::State& t_state) {
  const std::string name = SanitizeName(t_state.name());
  {
    std::ofstream ofs("construction-" + name + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
  }
  {
    std::ofstream ofs("construction-" + name + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
  }
}

void WriteSizesSidecar(const benchmark::State& t_state, const dret::SizeReport& sizes) {
  const std::string name = SanitizeName(t_state.name());
  dret::writeSizesJson("sizes-" + name + ".json", sizes);
}

//~~~~~~~  Query-mode result containers (vector of doc ids) ~~~~~~~

// Pure collecting sink. Under the set-semantics output contract, document
// listing returns the *set* of documents in unspecified order, and every index
// already delivers a distinct set from Search (RMQ dedups via its marked
// bitvector; the sampled-tree families sort+unique internally; Brute dedups in
// its own Search). So the benchmark neither sorts nor dedups in the timed loop
// — it just materialises the reported ids. (The correctness dump below applies
// its own sort+unique to get a canonical form for cross-index comparison.)
struct DocListResult {
  std::vector<std::size_t> docs;
  void operator()(std::size_t d) { docs.push_back(d); }
};


//~~~~~~~  Per-pattern timing (-per_pattern_dir) ~~~~~~~
//
// shannon is a Xen guest whose clocksource is `xen`, where steady_clock::now()
// costs ~1 us -- against 7-20 us for a whole query on page -- so it cannot
// bracket a single query. A TSC read costs ~11 ns there (constant_tsc), so the
// timed loop carries ~23 ns per pattern. Ticks are converted to ns with one
// factor calibrated against steady_clock at startup.
namespace per_pattern_timing {

inline std::uint64_t ReadTicks() {
#if defined(__x86_64__) || defined(__i386__)
  return __rdtsc();
#else
  return static_cast<std::uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count());
#endif
}

double CalibrateNsPerTick() {
#if defined(__x86_64__) || defined(__i386__)
  using Clock = std::chrono::steady_clock;
  const auto t0 = Clock::now();
  const auto c0 = ReadTicks();
  std::this_thread::sleep_for(std::chrono::milliseconds(250));
  const auto c1 = ReadTicks();
  const auto t1 = Clock::now();
  return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(c1 - c0);
#else
  return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::duration(1)).count();
#endif
}

// Shared by every benchmark of the run: the patterns and the text do not change
// between configurations, so occurrences are counted once, before timing starts.
struct Context {
  double ns_per_tick = 0;
  std::vector<std::size_t> occs;
  std::string patterns_path;
};

Context& Ctx() {
  static Context ctx;
  return ctx;
}

// Nearest-rank quantile of an ascending vector.
double Quantile(const std::vector<double>& t_sorted, double t_q) {
  if (t_sorted.empty()) return 0;
  const auto rank = static_cast<std::size_t>(t_q * static_cast<double>(t_sorted.size()));
  return t_sorted[std::min(rank, t_sorted.size() - 1)];
}

// Writes the sidecar and the quantile counters. Runs after the timed loop.
void Report(benchmark::State& t_state,
            const std::vector<std::string>& t_patterns,
            const std::vector<std::uint64_t>& t_ticks,
            const std::vector<std::size_t>& t_ndocs) {
  const auto& ctx = Ctx();
  const auto n = t_patterns.size();

  std::vector<double> ns(n);
  for (std::size_t i = 0; i < n; ++i) ns[i] = static_cast<double>(t_ticks[i]) * ctx.ns_per_tick;

  // Google Benchmark calls the benchmark function once per trial run (while it sizes
  // the iteration count) and once per repetition, so every call gets its own file;
  // the iteration count in the header tells trial runs from measured repetitions.
  static std::map<std::string, int> calls;
  const int call = ++calls[t_state.name()];
  const auto path = FLAGS_per_pattern_dir + "/" + SanitizeName(t_state.name()) + ".call" +
                    std::to_string(call) + ".tsv";
  std::ofstream os(path);
  if (!os) {
    std::cerr << "ERROR: cannot write per-pattern times to " << path << "\n";
  } else {
    os << "# benchmark\t" << t_state.name() << '\n'
       << "# call\t" << call << '\n'
       << "# patterns\t" << ctx.patterns_path << '\n'
       << "# pattern_idx counts the non-empty lines of the patterns file, from 0\n"
       << "# ns_per_tick\t" << std::setprecision(9) << ctx.ns_per_tick << '\n'
       << "# time_ns is the minimum over " << t_state.iterations() << " iteration(s)\n"
       << "pattern_idx\tlength\toccs\tndocs\ttime_ns\n"
       << std::fixed << std::setprecision(1);
    for (std::size_t i = 0; i < n; ++i)
      os << i << '\t' << t_patterns[i].size() << '\t' << ctx.occs[i] << '\t'
         << t_ndocs[i] << '\t' << ns[i] << '\n';
  }

  // A document cannot be reported without an occurrence in it.
  std::size_t ndocs_above_occs = 0;
  for (std::size_t i = 0; i < n; ++i) ndocs_above_occs += (t_ndocs[i] > ctx.occs[i]);
  if (ndocs_above_occs)
    std::cerr << "[VALIDATE] " << t_state.name() << ": " << ndocs_above_occs
              << " pattern(s) report more documents than occurrences\n";

  std::sort(ns.begin(), ns.end());
  t_state.counters["Time_x_Pattern_p50"] = Quantile(ns, 0.50) * 1e-9;
  t_state.counters["Time_x_Pattern_p90"] = Quantile(ns, 0.90) * 1e-9;
  t_state.counters["Time_x_Pattern_p99"] = Quantile(ns, 0.99) * 1e-9;
}

}  // namespace per_pattern_timing

//~~~~~~~  Counter helpers  ~~~~~~~

void SetupQueryCounters(benchmark::State& t_state) {
  t_state.counters["Size(bytes)"] = 0;
  t_state.counters["Bits_x_Symbol"] = 0;
  t_state.counters["Patterns"] = 0;
  t_state.counters["Time_x_Pattern"] = 0;
}

// The size and time counters of a query benchmark.
void QueryCounters(benchmark::State& t_state,
                   const dret::DocListIndex<>* idx,
                   std::size_t idx_size,
                   std::size_t t_n_patterns,
                   std::size_t t_seq_size) {
  SetupQueryCounters(t_state);

  // Headline size = the *resident* (in-memory) footprint from GetSizeReport,
  // which includes the rank/select/BP supports that serialize() rebuilds on
  // load (e.g. CompactBP's bp_support, CompactLOUDS's leaf_rank) and so are
  // absent from the serialized size. This keeps the cross-family Bits_x_Symbol
  // comparison apples-to-apples — families that rebuild supports on load would
  // otherwise be undercounted. The serialized (on-disk) size is still reported
  // as `serialized_bytes`. Brute baselines have no GetSizeReport (total == 0),
  // so they fall back to the serialized size.
  std::size_t resident = idx_size;
  if (idx) {
    const auto sizes = idx->GetSizeReport();
    appendCounters(t_state, sizes);
    const auto report_total = dret::totalBytes(sizes);
    t_state.counters["total_index_bytes"] = static_cast<double>(report_total);
    if (report_total > 0) resident = report_total;
  }

  t_state.counters["Size(bytes)"] = static_cast<double>(resident);
  t_state.counters["Bits_x_Symbol"] = resident * 8.0 / t_seq_size;
  t_state.counters["serialized_bytes"] = static_cast<double>(idx_size);
  t_state.counters["Patterns"] = t_n_patterns;
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_n_patterns,
      benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
}


//~~~~~~~  Query benchmark kernel ~~~~~~~

auto BM_Query = [](benchmark::State& t_state,
                   Factory<>* t_factory,
                   const Cell* t_cell,
                   const std::vector<std::string>* t_patterns,
                   std::size_t t_seq_size) {
  [[maybe_unused]] auto [idx, idx_size, frequencies] = t_factory->Get(*t_cell);

  // The per-pattern variant is a separate loop so the default path stays exactly
  // the plain one. It adds two TSC reads and two stores per pattern, into vectors
  // sized before the timed region; each pattern keeps its fastest iteration.
  const bool per_pattern = !FLAGS_per_pattern_dir.empty() && idx;
  std::vector<std::uint64_t> pattern_ticks;
  std::vector<std::size_t> pattern_ndocs;
  if (per_pattern) {
    pattern_ticks.assign(t_patterns->size(), std::numeric_limits<std::uint64_t>::max());
    pattern_ndocs.assign(t_patterns->size(), 0);
  }

#ifdef DRET_DOC_LIST_PROFILE
  // Phase profile of the sampled-tree indexes (dret/doc_list/search_profile.h), over the timed loop
  // only: read before the untimed passes below, which also search.
  dret::search_profile().reset();
#endif

  if (!per_pattern) {
    for (auto _ : t_state) {
      for (const auto& pattern : *t_patterns) {
        DocListResult result;
        idx->Search(pattern, std::ref(result));
        benchmark::DoNotOptimize(result);
      }
    }
  } else {
    for (auto _ : t_state) {
      for (std::size_t i = 0; i < t_patterns->size(); ++i) {
        const auto t0 = per_pattern_timing::ReadTicks();
        std::size_t n_docs;
        {
          DocListResult result;
          idx->Search((*t_patterns)[i], std::ref(result));
          benchmark::DoNotOptimize(result);
          n_docs = result.docs.size();
        }
        const auto t1 = per_pattern_timing::ReadTicks();
        pattern_ticks[i] = std::min(pattern_ticks[i], t1 - t0);
        pattern_ndocs[i] = n_docs;
      }
    }
  }

#ifdef DRET_DOC_LIST_PROFILE
  {
    const auto& pr = dret::search_profile();
    const double q = pr.n_queries ? static_cast<double>(pr.n_queries) : 1.0;
    t_state.counters["prof_count_ns_x_q"] = pr.ns_count / q;
    t_state.counters["prof_cover_ns_x_q"] = pr.ns_cover / q;
    t_state.counters["prof_expand_ns_x_q"] = pr.ns_expand / q;
    t_state.counters["prof_combine_ns_x_q"] = pr.ns_combine / q;
    t_state.counters["prof_nodes_x_q"] = pr.n_nodes / q;
    t_state.counters["prof_rawpos_x_q"] = pr.n_raw_positions / q;
    t_state.counters["prof_docs_x_q"] = pr.n_docs / q;
  }
#endif

  // Optional correctness dump: each index's per-pattern doc set, one line per
  // pattern, in canonical (sorted) order so identical answer sets compare equal
  // across families regardless of each index's native report order.
  //
  // We deliberately sort but do NOT deduplicate: every index must report each
  // document at most once, so erasing duplicates here would *hide* a buggy
  // index that double-reports a doc (it would still match the others). Keeping
  // duplicates makes such a bug surface as a file/md5 mismatch, and we also
  // count them for a direct, self-contained signal on stderr.
  if (!FLAGS_print_results_dir.empty() && idx) {
    const std::string fname = SanitizeName(t_state.name());
    std::ofstream os(FLAGS_print_results_dir + "/" + fname + ".txt");
    std::size_t dup_count = 0;
    std::size_t patterns_with_dups = 0;
    for (const auto& pattern : *t_patterns) {
      DocListResult result;
      idx->Search(pattern, std::ref(result));
      std::sort(result.docs.begin(), result.docs.end());
      std::size_t dups_here = 0;
      for (std::size_t i = 1; i < result.docs.size(); ++i)
        if (result.docs[i] == result.docs[i - 1]) ++dups_here;
      dup_count += dups_here;
      patterns_with_dups += (dups_here > 0);
      for (std::size_t i = 0; i < result.docs.size(); ++i)
        os << (i ? " " : "") << result.docs[i];
      os << '\n';
    }
    if (dup_count)
      std::cerr << "[VALIDATE] " << t_state.name() << ": reported " << dup_count
                << " duplicate doc-id(s) across " << patterns_with_dups
                << " pattern(s) (expected 0 — possible index bug)\n";
  }

  if (per_pattern) per_pattern_timing::Report(t_state, *t_patterns, pattern_ticks, pattern_ndocs);

  QueryCounters(t_state, idx.get(), idx_size, t_patterns->size(), t_seq_size);
};

// The documents of each pattern with their frequencies, for a frequency task. The answers dump
// writes, per pattern, "d:f" for each document in document order.
struct DocFreqResult {
  std::vector<std::pair<std::size_t, std::size_t>> freqs;
  void operator()(std::size_t d, std::size_t f) { freqs.emplace_back(d, f); }
};

auto BM_QueryFreq = [](benchmark::State& t_state,
                       Factory<>* t_factory,
                       const Cell* t_cell,
                       const std::vector<std::string>* t_patterns,
                       std::size_t t_seq_size) {
  auto [idx, idx_size, frequencies] = t_factory->Get(*t_cell);
  if (!frequencies) {
    t_state.SkipWithError("the index does not compute document frequencies");
    return;
  }
  for (auto _ : t_state) {
    for (const auto& pattern : *t_patterns) {
      DocFreqResult result;
      frequencies(pattern, std::ref(result));
      benchmark::DoNotOptimize(result);
    }
  }

  if (!FLAGS_print_results_dir.empty()) {
    std::ofstream os(FLAGS_print_results_dir + "/" + SanitizeName(t_state.name()) + ".txt");
    std::size_t dup_count = 0;
    for (const auto& pattern : *t_patterns) {
      DocFreqResult result;
      frequencies(pattern, std::ref(result));
      std::sort(result.freqs.begin(), result.freqs.end());
      for (std::size_t i = 0; i < result.freqs.size(); ++i) {
        if (i && result.freqs[i].first == result.freqs[i - 1].first) ++dup_count;
        os << (i ? " " : "") << result.freqs[i].first << ':' << result.freqs[i].second;
      }
      os << '\n';
    }
    if (dup_count)
      std::cerr << "[VALIDATE] " << t_state.name() << ": reported " << dup_count
                << " document(s) twice (expected 0 — possible index bug)\n";
  }

  QueryCounters(t_state, idx.get(), idx_size, t_patterns->size(), t_seq_size);
};

//~~~~~~~  Construct benchmark kernels ~~~~~~~
//
// Construct kernels build a fresh, owning index and time the construct() call.
// The benchmark instance owns a GenericStorage so the typed index isn't tied
// to the factory's external-storage cache.

void SetupConstructCounters(benchmark::State& t_state, const dret::Config& t_config,
                             std::uint32_t bs, float sf) {
  using namespace sri::conf;
  t_state.counters["n"] = ReadVectorLength(sdsl::cache_file_name(t_config.keys[kBWT][kBase], t_config));
  if (bs > 0) t_state.counters["bs"] = bs;
  if (sf > 0) t_state.counters["sf"] = sf;
}


// One kernel for every cell: the index is made once, in a storage the benchmark owns; each iteration
// times its construct() (with --rebuild, after deleting the cell's own cache files, untimed); then it
// is loaded once for the size report.
void BM_Construct(benchmark::State& t_state,
                  const Cell* t_cell,
                  dret::Config* t_config,
                  std::function<void()> t_rebuild_hook,
                  bool t_memory_trace) {
  sri::GenericStorage storage;
  auto cell_index = t_cell->make(std::ref(storage));
  std::int64_t first_iter_ns = -1;
  for (auto _ : t_state) {
    if (t_rebuild_hook) {
      t_state.PauseTiming();
      t_rebuild_hook();
      t_state.ResumeTiming();
    }
    const auto t0 = std::chrono::steady_clock::now();
    sdsl::memory_monitor::start();
    cell_index.construct(*t_config);
    sdsl::memory_monitor::stop();
    const auto t1 = std::chrono::steady_clock::now();
    if (first_iter_ns < 0) first_iter_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
  }
  WarnIfWarm(t_state, first_iter_ns, static_cast<bool>(t_rebuild_hook));
  if (t_memory_trace) WriteMemoryTrace(t_state);
  SetupConstructCounters(t_state, *t_config, t_cell->block_size, t_cell->storing_factor);
  if (t_cell->sampling > 0) t_state.counters["s"] = static_cast<double>(t_cell->sampling);
  t_state.counters["first_construct_ns"] = static_cast<double>(first_iter_ns);

  cell_index.load(*t_config);
  const auto sizes = cell_index.index->GetSizeReport();
  appendCounters(t_state, sizes);
  t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  if (t_memory_trace) WriteSizesSidecar(t_state, sizes);
}

// The cells of the spec, in its order. They outlive the benchmarks, which point to them.
std::vector<Cell> SpecCells(const bench::spec::Spec& t_spec) {
  std::vector<Cell> cells;
  for (const auto& block : t_spec.sweep) {
    auto block_cells = std::visit([](const auto& t_sweep) { return bench::factories::Cells(t_sweep); }, block);
    for (auto& c : block_cells) cells.push_back(std::move(c));
  }
  return cells;
}

}  // namespace

int main(int argc, char** argv) {
  // Stamp the vocabulary these labels are written in. experiments/analyze.py
  // refuses a file without it, and migrate_result_names.py skips a file that
  // has it -- which is what stops the migrator renaming an already-correct
  // fresh run (its rewrite is not idempotent: SADA would become SADA-L).
  benchmark::AddCustomContext("core_naming", "2026-11");

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

  // Carry the spec's RePair block into the config the construct code sees. It is
  // deliberately not read from the environment: <MB> selects between two
  // different grammars, so the value has to be part of the archived spec rather
  // than ambient process state.
  {
    const auto& r = spec.repair;
    config.repair.variant = r.variant == "bal32"   ? dret::repair::Variant::kBal32
                            : r.variant == "bal64" ? dret::repair::Variant::kBal64
                                                   : dret::repair::Variant::kAuto;
    config.repair.mb = r.mb;
    config.repair.mb_max = r.mb_max;
    if (spec.mode == bench::spec::Mode::Construct &&
        (r.variant != "auto" || r.mb.has_value() || r.mb_max.has_value())) {
      std::cerr << "RePair: spec override in effect -- variant=" << r.variant
                << " mb=" << (r.mb ? std::to_string(*r.mb) : "derived")
                << " mb_max=" << (r.mb_max ? std::to_string(*r.mb_max) : "none") << std::endl;
    }
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

    if (!FLAGS_per_pattern_dir.empty()) {
      std::filesystem::create_directories(FLAGS_per_pattern_dir);
      auto& ctx = per_pattern_timing::Ctx();
      ctx.patterns_path = spec.dataset.patterns;
      ctx.ns_per_tick = per_pattern_timing::CalibrateNsPerTick();
      ctx.occs.reserve(patterns.size());
      for (const auto& pattern : patterns) ctx.occs.push_back(factory->CountOccurrences(pattern));
      std::cerr << "per-pattern timing: " << patterns.size() << " patterns, ns/tick = "
                << ctx.ns_per_tick << ", writing to " << FLAGS_per_pattern_dir << std::endl;
    }
  }

  // One benchmark per cell, the same names in both modes.
  const auto cells = SpecCells(spec);
  for (const auto& cell : cells) {
    if (spec.mode == bench::spec::Mode::Query && spec.workload.task == bench::spec::Task::Frequency) {
      benchmark::RegisterBenchmark(cell.name, BM_QueryFreq, factory.get(), &cell, &patterns, seq_size);
    } else if (spec.mode == bench::spec::Mode::Query) {
      benchmark::RegisterBenchmark(cell.name, BM_Query, factory.get(), &cell, &patterns, seq_size);
    } else {
      // In construct mode the cache directory is the working directory (the dret::Config above) and
      // the cache files end in _<dataset name>.sdsl.
      auto hook = cache_clean::MakeHook(
          spec.workload.rebuild, std::filesystem::current_path(), spec.dataset.name, cell.rebuild_prefixes);
      benchmark::RegisterBenchmark(cell.name, BM_Construct, &cell, &config, hook, spec.workload.memory_trace);
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  return 0;
}
