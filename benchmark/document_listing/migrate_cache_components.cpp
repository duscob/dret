//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/4/26.
//
// Converts a cache directory written before 2026-10 into cache components (see
// include/dret/cache_components.h), without rebuilding anything:
//
//   - files whose format did not change get their new name as a hard link (the
//     DA grammars, GCDA-nolists' plain grammar, the node document lists);
//   - objects now stored as components are loaded and stored as them (every
//     GCDA representation, GCDA-differential, GCDA-nolists' differential
//     grammar, the PDL cores), then loaded back from the components and checked
//     against the original;
//   - build intermediates no index reads any more (the differential base-grammar
//     caches, the raw LightSLP<>) are obsolete.
//
// Old files are left in place unless --delete_old, which removes the old names
// of everything converted and the obsolete files. Files of other kinds (text,
// SA, DA, the r-index, the RMQ cores) keep their names and are not touched.
//
//   migrate_cache_components --dir=<cache dir> --id=<collection id> [--dry_run] [--delete_old]
//

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <regex>
#include <string>
#include <vector>

#include <gflags/gflags.h>

#include "dret/cache_components.h"
#include "dret/config.h"

#include "factories/dgcda.h"
#include "factories/gcda.h"
#include "factories/pdl.h"
#include "factories/slp_ns.h"

DEFINE_string(dir, "", "Cache directory of one collection (MANDATORY).");
DEFINE_string(id, "", "Collection id: the suffix of its cache files, <key>_<type>_<id>.sdsl (MANDATORY).");
DEFINE_bool(dry_run, false, "Report what would be done; write nothing.");
DEFINE_bool(delete_old, false, "Remove the old names of converted files, and the obsolete files.");

namespace {

namespace fs = std::filesystem;
using dret::Component;
using dret::TypeHash;

struct OldFile {
  fs::path path;
  std::string key;
  std::string type;
};

struct Stats {
  std::size_t linked = 0, split = 0, obsolete = 0, unknown = 0, skipped = 0, failed = 0;
  std::vector<fs::path> to_delete;
};

const dret::JSON& Keys() { return dret::kDefaultKeys.keys; }

sdsl::cache_config CacheConfig() {
  sdsl::cache_config config(false, FLAGS_dir, FLAGS_id);
  return config;
}

std::optional<OldFile> Parse(const fs::path& t_path) {
  const auto name = t_path.filename().string();
  const auto suffix = "_" + FLAGS_id + ".sdsl";
  if (name.size() <= suffix.size() || !name.ends_with(suffix)) return std::nullopt;
  const auto stem = name.substr(0, name.size() - suffix.size());
  const auto cut = stem.rfind('_');
  if (cut == std::string::npos) return std::nullopt;
  const auto type = stem.substr(cut + 1);
  if (type.empty() || type.find_first_not_of("0123456789") != std::string::npos) return std::nullopt;
  return OldFile{t_path, stem.substr(0, cut), type};
}

void Link(const OldFile& t_old, const std::string& t_new_key, Stats& t_stats) {
  const auto target = dret::ComponentFile({t_new_key, t_old.type}, CacheConfig());
  if (fs::exists(target)) {
    ++t_stats.skipped;
  } else {
    std::cout << "link  " << t_old.path.filename().string() << " -> " << fs::path(target).filename().string() << "\n";
    if (!FLAGS_dry_run) fs::create_hard_link(t_old.path, target);
    ++t_stats.linked;
  }
  t_stats.to_delete.push_back(t_old.path);
}

// Load the old object with t_load, store it as components, load it back from
// them and compare. A type with operator== must compare equal; every type must
// come back at the same size.
template <typename T, typename TLoad, typename TComponents>
bool Split(const OldFile& t_old, TLoad&& t_load, TComponents&& t_components, Stats& t_stats) {
  T object;
  {
    std::ifstream in(t_old.path, std::ios::binary);
    t_load(object, in);
    if (!in) throw std::runtime_error("cannot read " + t_old.path.string());
  }
  const auto components = t_components(object);
  std::cout << "split " << t_old.path.filename().string() << " ->";
  for (const auto& c : components) std::cout << " " << c.key;
  std::cout << "\n";
  t_stats.to_delete.push_back(t_old.path);
  if (FLAGS_dry_run) return ++t_stats.split, true;

  dret::StoreComponents(object, components, CacheConfig());
  T back;
  if (!dret::LoadComponents(back, components, CacheConfig()) ||
      dret::SerializedSize(back) != dret::SerializedSize(object)) {
    std::cerr << "FAILED to read back " << t_old.path << "\n";
    t_stats.to_delete.pop_back();
    return ++t_stats.failed, false;
  }
  if constexpr (requires { back == object; }) {
    if (!(back == object)) {
      std::cerr << "FAILED: " << t_old.path << " reads back as another object\n";
      t_stats.to_delete.pop_back();
      return ++t_stats.failed, false;
    }
  }
  ++t_stats.split;
  return true;
}

// Dispatch t_old over the candidate types T...: the one whose type hash it carries.
template <typename... T, typename TFn>
bool ForType(const OldFile& t_old, TFn&& t_fn) {
  return ((t_old.type == TypeHash<T>() ? (t_fn.template operator()<T>(), true) : false) || ...);
}

template <typename T>
auto LoadCurrent() {
  return [](T& t, std::istream& in) { t.load(in); };
}

template <typename T>
auto LoadPre202610() {
  return [](T& t, std::istream& in) { t.loadPre202610(in); };
}

}  // namespace

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_dir.empty() || FLAGS_id.empty()) {
    std::cerr << "usage: migrate_cache_components --dir=<cache dir> --id=<collection id> [--dry_run] [--delete_old]\n";
    return 1;
  }

  namespace gcda = bench::factories::gcda;
  namespace dgcda = bench::factories::dgcda;
  namespace slp_ns = bench::factories::slp_ns;
  using namespace dret::conf;
  using PDLCore = dret::pdl::DocListIdxPDLPlain<dret::GenericStorage>::TCore;
  using CorePlain = typename dret::pdl::WithCodec<PDLCore, dret::pdl::PlainCodec<>>::type;
  using CoreRP = typename dret::pdl::WithCodec<PDLCore, dret::pdl::RPCodec<>>::type;
  using CoreRPLengths = typename dret::pdl::WithCodec<PDLCore, dret::pdl::RPCodecWithLengths>::type;
  using CoreBC = typename dret::pdl::WithCodec<PDLCore, dret::pdl::BCCodec<>>::type;

  const std::regex cell_re(R"((\d+)-([0-9.]+)_(gcda_slp|gcda_docs|dgcda_slp|dgcda_docs))");
  const std::regex pdl_re(R"((\d+)-([0-9.]+)_pdl_(plain|rp|bc)_(\d)_pdl_(plain|rp|bc)_sets)");
  const std::regex dslp_re(R"((?:bs(\d+)_)?dslp_ns)");
  const std::regex obsolete_re(R"(dgcda_slp_grammar|(?:bs\d+_)?dslp_ns_grammar|bs\d+_slp_ns(?:_grammar)?)");

  Stats stats;
  std::vector<OldFile> files;
  for (const auto& e : fs::directory_iterator(FLAGS_dir))
    if (auto f = Parse(e.path())) files.push_back(*f);
  std::sort(files.begin(), files.end(), [](const auto& a, const auto& b) { return a.path < b.path; });

  for (const auto& f : files) {
    std::smatch m;
    try {
      if (f.key == "gcda_slp_grammar") {
        Link(f, dret::KeyName(Keys(), kDaGrammar), stats);
      } else if (f.key == "gcda_slp_compact_seq") {
        Link(f, dret::KeyName(Keys(), kDaGrammarSeq), stats);
      } else if (f.key == "dgcda_slp_cnf" || f.key == "slp_ns") {
        // The CNF grammar: DGCDA's raw copy, and GCDA-nolists' plain grammar in
        // its containers. Neither format changed.
        Link(f, dret::KeyName(Keys(), kDaCnfGrammar), stats);
      } else if (std::regex_match(f.key, m, obsolete_re)) {
        std::cout << "obsolete " << f.path.filename().string() << "\n";
        stats.to_delete.push_back(f.path);
        ++stats.obsolete;
      } else if (std::regex_match(f.key, m, cell_re)) {
        const dret::SampledTreeCell cell{static_cast<uint32_t>(std::stoul(m[1])), std::stof(m[2])};
        const auto what = m[3].str();
        bool known = false;
        if (what == "gcda_docs" || what == "dgcda_docs") {
          known = ForType<grammar::Chunks<>, gcda::Sets_Plain>(f, [&]<typename T>() {
                    Link(f, dret::CellKey(Keys(), kDaNodeDocListsPlain, cell), stats);
                  }) ||
                  ForType<grammar::GCChunks<grammar::SLP<>>, gcda::Sets_RP>(f, [&]<typename T>() {
                    Link(f, dret::CellKey(Keys(), kDaNodeDocListsRP, cell), stats);
                  });
        } else if (what == "gcda_slp") {
          auto split = [&]<typename T>() {
            Split<T>(f, LoadCurrent<T>(),
                     [&](const T& t) { return dret::CacheComponents(t, Keys(), cell); }, stats);
          };
          known = ForType<grammar::CombinedSLP<>, gcda::SLP_Light, gcda::SLP_Combined, gcda::SLP_CompactBP,
                          gcda::SLP_CompactLOUDS>(f, split);
          if (!known && f.type == TypeHash<grammar::LightSLP<>>()) {
            // The raw LightSLP<> was a build intermediate; nothing reads it now.
            std::cout << "obsolete " << f.path.filename().string() << "\n";
            stats.to_delete.push_back(f.path);
            ++stats.obsolete;
            known = true;
          }
        } else {  // dgcda_slp
          known = ForType<dgcda::SLP_Default, dgcda::SLP_OTF, dgcda::SLP_CRL, dgcda::SLP_EV, dgcda::SLP_DV,
                          dgcda::SLP_VV>(f, [&]<typename T>() {
            Split<T>(f, LoadPre202610<T>(),
                     [&](const T& t) { return dret::CacheComponents(t, Keys(), cell); }, stats);
          });
        }
        if (!known) {
          std::cout << "unknown type, left alone: " << f.path.filename().string() << "\n";
          ++stats.unknown;
        }
      } else if (std::regex_match(f.key, m, dslp_re)) {
        const uint32_t spacing = m[1].matched ? static_cast<uint32_t>(std::stoul(m[1])) : dret::kDiffBlockSize;
        const bool known = ForType<slp_ns::BareSLP_Diff, slp_ns::BareSLP_DiffEV, slp_ns::BareSLP_DiffDV,
                                   slp_ns::BareSLP_DiffVV>(f, [&]<typename T>() {
          Split<T>(f, LoadPre202610<T>(),
                   [&](const T& t) { return dret::CacheComponents(t, Keys(), spacing); }, stats);
        });
        if (!known) {
          std::cout << "unknown type, left alone: " << f.path.filename().string() << "\n";
          ++stats.unknown;
        }
      } else if (std::regex_match(f.key, m, pdl_re)) {
        const auto b = static_cast<uint32_t>(std::stoul(m[1]));
        const float sf = std::stof(m[2]);
        const auto policy = static_cast<dret::pdl::StoragePolicy>(std::stoi(m[4]));
        auto components = [&](const auto& core) { return dret::pdl::CacheComponents(core, Keys(), b, sf, policy); };
        bool known = ForType<CorePlain, CoreRP, CoreBC>(f, [&]<typename T>() {
          Split<T>(f, LoadPre202610<T>(), components, stats);
        });
        if (!known && f.type == TypeHash<CoreRPLengths>()) {
          // Re-Pair lists cached with the per-rule length array they never read:
          // convert, as construct() used to.
          Split<CoreRP>(f, [](CoreRP& t, std::istream& in) {
            CoreRPLengths legacy;
            legacy.loadPre202610(in);
            t = CoreRP(legacy);
          }, components, stats);
          known = true;
        }
        if (!known) {
          std::cout << "unknown type, left alone: " << f.path.filename().string() << "\n";
          ++stats.unknown;
        }
      }
    } catch (const std::exception& e) {
      std::cerr << "FAILED " << f.path << ": " << e.what() << "\n";
      ++stats.failed;
    }
  }

  std::cout << "linked " << stats.linked << ", split " << stats.split << ", obsolete " << stats.obsolete
            << ", already there " << stats.skipped << ", unknown " << stats.unknown << ", failed " << stats.failed
            << "\n";
  if (FLAGS_delete_old && !FLAGS_dry_run) {
    if (stats.failed) {
      std::cerr << "not deleting old files: " << stats.failed << " failures\n";
      return 2;
    }
    for (const auto& p : stats.to_delete) fs::remove(p);
    std::cout << "deleted " << stats.to_delete.size() << " old files\n";
  }
  return stats.failed ? 2 : 0;
}
