//
// JSON sweep-spec types + parser for bm_doc_list.
//
// A spec describes one benchmark run: the dataset, the workload settings,
// the output destinations, the mode (query vs construct), and a list of
// per-family sweep blocks. Each block enumerates the axes a single family
// cares about; the driver expands the Cartesian product per block.
//
// Schema (informal):
//
//   {
//     "spec_version": 2,    (the RMQ core names; see ParseRmqCore)
//     "mode":     "query" | "construct",
//     "dataset":  { "dir": "...", "name": "...", "patterns": "...",
//                   "data_width": 8, "sa_algo": "SDSL_LIBDIVSUFSORT",
//                   "doc_delim": 0 },
//     "workload": { "reps": 10, "min_time": 0.0,
//                   "report_stats": false, "print_result": false,
//                   "task": "listing" | "frequency" },
//     "repair":   { "variant": "auto", "mb": null, "mb_max": null },
//     "sweep": [
//       { "family": "gcda",  "tslp": ["Light"], "block_size": [512],
//                            "storing_factor": [4] },
//       { "family": "dgcda", "tslp": ["Default","OTF"], "block_size": [512],
//                            "storing_factor": [4] },
//       { "family": "rmq",   "core": ["sada","ilcp","cilcp-v"],   (needs "spec_version": 2)
//                            "get_doc": ["da","slp"], "gcda_slp": ["Light"],
//                            "block_size": [512], "storing_factor": [4] },
//       { "family": "pdl",   "codec": ["Plain"], "get_doc": ["DA"],
//                            "policy": ["OccurrenceWeighted"],
//                            "block_size": [512], "storing_factor": [4] },
//       { "family": "slp_ns","tslp": ["IV"] },
//       { "family": "brute", "kind": ["r-index"], "sampling_size": [16] }
//       { "family": "brute", "kind": ["da", "rlcsa"] }   (brute force over the DA / RLCSA backend)
//     ]
//   }
//
// Unknown family keys throw; missing optional keys default to the
// corresponding axis enum's TryParse default value.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#include <nlohmann/json.hpp>

#include "dret/config.h"

#include "axes.h"
#include "enum_traits.h"

namespace bench::spec {

enum class Mode { Query, Construct };

// What a query returns: the documents of the pattern, or each with its frequency (query mode only).
enum class Task { Listing, Frequency };

struct Dataset {
  std::string dir;
  std::string name;
  std::string patterns;  // Query mode only.
  std::uint8_t data_width = 8;
  std::string sa_algo = "SDSL_LIBDIVSUFSORT";
  int doc_delim = 0;
};

// How to run irepair over this collection's document array. Every field is
// optional and the defaults derive from the array's size, so a spec that omits
// the block builds exactly what it always did.
//
// It belongs in the spec rather than in an environment variable because <MB>
// decides *which grammar gets built*: below a size-dependent threshold irepair
// switches to a memory-lean pass whose output differs. The spec is archived with
// the run, so a value set here leaves a record; an environment variable would
// not. Construct mode only.
struct RePair {
  // "auto" (default) | "bal32" | "bal64". Both are RePair's balanced build; the
  // number is the width of a sequence position, and only bal64 is correct past
  // 2^31 elements. See include/dret/repair.h.
  std::string variant = "auto";
  // Unset: derive the value that keeps irepair on its fast path.
  std::optional<std::uint64_t> mb;
  // Unset: no cap. Setting this below the derived value trades the grammar for
  // memory, and the construct logs say so when it happens.
  std::optional<std::uint64_t> mb_max;
};

struct Workload {
  int reps = 10;
  double min_time = 0.0;
  bool report_stats = false;
  bool print_result = false;
  // Construct mode only: delete the cell's variant-specific cache files
  // before each timed construct() iteration. The shared cache artefacts
  // (Text, SA, DocEnds, DA, LCP) are kept. Supported for GCDA, DGCDA,
  // SLP-NS, PDL — see bm_doc_list.cpp for the per-family glob patterns.
  // Default false (warm-cache fast path).
  bool rebuild = false;
  // Construct mode only: write construction-<bench>.{html,json} memory-monitor
  // traces and sizes-<bench>.json size sidecars next to the working directory
  // for every benchmark cell.
  bool memory_trace = false;
  // Pattern coding for query mode: "PLAIN" | "BASE64". Parsed but not yet read:
  // bm_doc_list reads patterns as PLAIN lines.
  std::string pattern_code = "PLAIN";
  int pattern_delim = '\n';
  Task task = Task::Listing;
};

struct GCDASweep {
  std::vector<axes::GCDASLPVariant> tslp{axes::GCDASLPVariant::Light};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
  // Document-set codec: "rp" (default) or "plain".
  std::vector<axes::GCDASetsCodec> sets{axes::GCDASetsCodec::RP};
};

struct DGCDASweep {
  std::vector<axes::DGCDASLPVariant> tslp{axes::DGCDASLPVariant::Default};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
  std::vector<axes::GCDASetsCodec> sets{axes::GCDASetsCodec::RP};
};

struct SLPNSSweep {
  std::vector<axes::BareSLPVariant> tslp{axes::BareSLPVariant::IV};
  // Sample spacing of the differential variants (plain ones have no sampling),
  // and block size of the sampled ones. Default: the spacing GCDA-nolists has
  // always used.
  std::vector<std::uint32_t> block_size{512};
  // Sampled GCDA grammars without their document lists, one per
  // (block_size, storing_factor) cell. None by default.
  std::vector<axes::GCDASLPVariant> gcda_slp{};
  std::vector<float> storing_factor{32.0f};
};

struct RMQSweep {
  std::vector<axes::CoreKind> core{axes::CoreKind::SADA,
                                   axes::CoreKind::ILCP,
                                   axes::CoreKind::CILCP,
                                   axes::CoreKind::SADA_V,
                                   axes::CoreKind::ILCP_V,
                                   axes::CoreKind::CILCP_V};
  std::vector<axes::GetDocEnum> get_doc{axes::GetDocEnum::DA};
  std::vector<axes::GCDASLPVariant> gcda_slp{axes::GCDASLPVariant::Light};
  std::vector<axes::BareSLPVariant> bare_slp{axes::BareSLPVariant::IV};
  // DGCDA variant axis for get_doc=dslp (the differential get-doc backend).
  std::vector<axes::DGCDASLPVariant> dgcda_slp{axes::DGCDASLPVariant::Default};
  // TRunValues axis for the value-rule ILCP-V / CILCP-V cores. Default DV
  // matches IlcpLikeValueCore's default; other variants in {iv, dv, vv} fan
  // the -V cores out across the run-values container.
  std::vector<axes::RunValuesVariant> run_values{axes::RunValuesVariant::DV};
  // TPrevDoc axis for the value-rule SADA-V core. Default IV matches its default;
  // other variants in {iv, dv, vv}.
  std::vector<axes::PrevDocVariant> prev_doc{axes::PrevDocVariant::IV};
  // SA-Phi sampling rate axis for get_doc=sa_phi_sr. Ignored for sa_phi_r
  // (r-index has no sampling knob) and for other get-docs.
  std::vector<std::size_t> sa_sampling{8};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
};

struct PDLSweep {
  std::vector<axes::PDLVariant> codec;  // No default — PDL disabled if list empty.
  std::vector<axes::GetDocEnum> get_doc{axes::GetDocEnum::DA};
  std::vector<axes::PDLStoragePolicy> policy{axes::PDLStoragePolicy::OccurrenceWeighted};
  // GCDA-family backend sub-axes, consulted per get_doc kind: gcda_slp for
  // get_doc=slp, bare_slp for get_doc=slp_ns, dgcda_slp for get_doc=dslp.
  // Defaults reproduce the pre-existing single-backend PDL behaviour.
  std::vector<axes::GCDASLPVariant> gcda_slp{axes::GCDASLPVariant::Light};
  std::vector<axes::BareSLPVariant> bare_slp{axes::BareSLPVariant::IV};
  std::vector<axes::DGCDASLPVariant> dgcda_slp{axes::DGCDASLPVariant::Default};
  // SA-Phi sampling rate axis for get_doc=sa_phi_sr (same shape as RMQSweep).
  std::vector<std::size_t> sa_sampling{8};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
  // Block size of the GCDA backend for get_doc=slp (gcda_slp light; its storing
  // factor is fixed at 32) and sample spacing for get_doc=slp_ns differential
  // variants. Default: the value the light backend used to be pinned at.
  std::vector<std::uint32_t> backend_block_size{1024};
};

// GCDA over a document-array backend (dret::gcda::DocListIdxGCDABackend): the
// sampled tree and node lists of GCDA, the ends of a range expanded by get_doc.
struct GCDABackendSweep {
  std::vector<axes::GetDocEnum> get_doc{axes::GetDocEnum::SAPhiR};
  std::vector<axes::GCDASetsCodec> sets{axes::GCDASetsCodec::RP};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
};

struct BruteSweep {
  enum class Kind { RIndex, SrIndex, DA, RLCSA, Count };
  std::vector<Kind> kind{Kind::RIndex};
  std::vector<std::size_t> sampling_size{8};  // Only used for sr-index.
};

using FamilySweep = std::variant<GCDASweep, DGCDASweep, SLPNSSweep,
                                  RMQSweep, PDLSweep, BruteSweep, GCDABackendSweep>;

struct Spec {
  // 0 when the spec does not say (see ParseRmqCore).
  int spec_version = 0;
  Mode mode = Mode::Query;
  Dataset dataset;
  Workload workload;
  RePair repair;
  std::vector<FamilySweep> sweep;
};

//~~~~~~~  Parser  ~~~~~~~

namespace detail {

inline Mode ParseMode(const std::string& s) {
  if (s == "query")     return Mode::Query;
  if (s == "construct") return Mode::Construct;
  throw std::invalid_argument("spec.mode: must be 'query' or 'construct', got '" + s + "'");
}

// Core vocabulary, versioned by spec.spec_version, because the names were
// swapped twice:
//   2026-09 (version 1): "sada" / "ilcp" / "cilcp" were the value-rule cores and
//     "sada-l" / "ilcp-l" / "cilcp-l" the marking ones; before that, the value
//     cores were the "-s" spellings and the plain names the marking cores.
//   2026-11 (version 2, current): the plain names are the published algorithms,
//     which stop by marking (Sadakane 2007; Gagie et al. 2017; CMR20), and
//     "sada-v" / "ilcp-v" / "cilcp-v" are dret's value-rule variants.
// A spec without spec_version may use only the unambiguous names ("-l" for the
// marking cores, "-v" for the value cores): a plain name there could mean either
// core, and reading it the wrong way would run the other core with plausible
// results. The "-s" spellings are rejected in every version.
inline constexpr int kSpecVersion = 2;

inline axes::CoreKind ParseRmqCore(const std::string& s, int t_version) {
  using axes::CoreKind;
  const bool plain = s == "sada" || s == "ilcp" || s == "cilcp";
  const auto base = [&s](std::size_t t_suffix) { return s.substr(0, s.size() - t_suffix); };
  auto marking = [](const std::string& b) {
    return b == "sada" ? CoreKind::SADA : b == "ilcp" ? CoreKind::ILCP : CoreKind::CILCP;
  };
  auto value = [](const std::string& b) {
    return b == "sada" ? CoreKind::SADA_V : b == "ilcp" ? CoreKind::ILCP_V : CoreKind::CILCP_V;
  };
  auto suffixed = [&s](const char* t_suffix) {
    for (const char* b : {"sada", "ilcp", "cilcp"})
      if (s == std::string(b) + t_suffix) return true;
    return false;
  };

  if (suffixed("-s")) {
    throw std::invalid_argument("rmq.core: '" + s + "' is a pre-2026-09 name; use '" + base(2) +
                                "-v' (the value-rule core) with spec_version " +
                                std::to_string(kSpecVersion));
  }
  if (plain) {
    if (t_version == 2) return marking(s);
    if (t_version == 1) return value(s);
    throw std::invalid_argument("rmq.core: '" + s + "' means the marking core under spec_version 2 and the "
                                "value-rule core under spec_version 1; set spec_version");
  }
  if (suffixed("-l")) {
    if (t_version == 2)
      throw std::invalid_argument("rmq.core: '" + s + "' is a 2026-09 name; under spec_version 2 the "
                                  "marking core is '" + base(2) + "'");
    return marking(base(2));
  }
  if (suffixed("-v")) {
    if (t_version == 1)
      throw std::invalid_argument("rmq.core: '" + s + "' does not exist under spec_version 1, where the "
                                  "value-rule core is '" + base(2) + "'");
    return value(base(2));
  }
  throw std::invalid_argument("rmq.core: must be 'sada' / 'ilcp' / 'cilcp' (marking) or 'sada-v' / "
                              "'ilcp-v' / 'cilcp-v' (value rule), got '" + s + "'");
}

inline BruteSweep::Kind ParseBruteKind(const std::string& s) {
  if (s == "r-index" || s == "rindex")    return BruteSweep::Kind::RIndex;
  if (s == "sr-index" || s == "srindex")  return BruteSweep::Kind::SrIndex;
  if (s == "da")                          return BruteSweep::Kind::DA;
  if (s == "rlcsa")                       return BruteSweep::Kind::RLCSA;
  if (s == "count")                       return BruteSweep::Kind::Count;
  throw std::invalid_argument("brute.kind: must be 'r-index', 'sr-index', 'da', 'rlcsa' or 'count', got '" + s + "'");
}

// Read a list-of-strings field and parse each via EnumTraits<E>::TryParse.
// Returns the default-initialised vector (the struct's in-class init) if the
// field is missing, so the spec author can omit irrelevant axes.
template <typename E>
std::vector<E> ParseEnumList(const nlohmann::json& j, const char* key,
                              const std::vector<E>& fallback) {
  if (!j.contains(key)) return fallback;
  std::vector<E> out;
  for (const auto& v : j.at(key)) {
    auto parsed = axes::EnumTraits<E>::TryParse(v.get<std::string>());
    if (!parsed) {
      throw std::invalid_argument(std::string(key) + ": unknown value '"
                                  + v.get<std::string>() + "'");
    }
    out.push_back(*parsed);
  }
  return out;
}

template <typename T, typename Parse>
std::vector<T> ParseListWith(const nlohmann::json& j, const char* key,
                              const std::vector<T>& fallback, Parse p) {
  if (!j.contains(key)) return fallback;
  std::vector<T> out;
  for (const auto& v : j.at(key)) out.push_back(p(v.get<std::string>()));
  return out;
}

template <typename T>
std::vector<T> ParseList(const nlohmann::json& j, const char* key,
                          const std::vector<T>& fallback) {
  if (!j.contains(key)) return fallback;
  return j.at(key).get<std::vector<T>>();
}

// Storing factors: numbers, or "inf" for a rule that stores no internal node
// (GCDA: no node's list beats its children's; PDL: no node's occurrences exceed
// infinity times its documents). JSON has no infinity, hence the string.
inline std::vector<float> ParseStoringFactors(const nlohmann::json& j, const std::vector<float>& fallback) {
  if (!j.contains("storing_factor")) return fallback;
  std::vector<float> out;
  for (const auto& v : j.at("storing_factor")) {
    if (v.is_string()) {
      if (v.get<std::string>() != "inf")
        throw std::invalid_argument("storing_factor: a number or \"inf\", not \"" + v.get<std::string>() + "\"");
      out.push_back(dret::kInfiniteStoringFactor);
    } else {
      out.push_back(v.get<float>());
    }
  }
  return out;
}

inline GCDASweep ParseGCDA(const nlohmann::json& j) {
  GCDASweep s;
  s.tslp = ParseEnumList<axes::GCDASLPVariant>(j, "tslp", s.tslp);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseStoringFactors(j, s.storing_factor);
  s.sets = ParseEnumList<axes::GCDASetsCodec>(j, "sets", s.sets);
  return s;
}

inline DGCDASweep ParseDGCDA(const nlohmann::json& j) {
  DGCDASweep s;
  s.tslp = ParseEnumList<axes::DGCDASLPVariant>(j, "tslp", s.tslp);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseStoringFactors(j, s.storing_factor);
  s.sets = ParseEnumList<axes::GCDASetsCodec>(j, "sets", s.sets);
  return s;
}

inline SLPNSSweep ParseSLPNS(const nlohmann::json& j) {
  SLPNSSweep s;
  s.tslp = ParseEnumList<axes::BareSLPVariant>(j, "tslp", s.tslp);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.gcda_slp = ParseEnumList<axes::GCDASLPVariant>(j, "gcda_slp", s.gcda_slp);
  s.storing_factor = ParseStoringFactors(j, s.storing_factor);
  return s;
}

inline RMQSweep ParseRMQ(const nlohmann::json& j, int t_version) {
  RMQSweep s;
  s.core = ParseListWith<axes::CoreKind>(
      j, "core", s.core, [t_version](const std::string& c) { return ParseRmqCore(c, t_version); });
  s.get_doc = ParseEnumList<axes::GetDocEnum>(j, "get_doc", s.get_doc);
  s.gcda_slp = ParseEnumList<axes::GCDASLPVariant>(j, "gcda_slp", s.gcda_slp);
  s.bare_slp = ParseEnumList<axes::BareSLPVariant>(j, "bare_slp", s.bare_slp);
  s.dgcda_slp = ParseEnumList<axes::DGCDASLPVariant>(j, "dgcda_slp", s.dgcda_slp);
  s.run_values = ParseEnumList<axes::RunValuesVariant>(j, "run_values", s.run_values);
  s.prev_doc = ParseEnumList<axes::PrevDocVariant>(j, "prev_doc", s.prev_doc);
  s.sa_sampling = ParseList<std::size_t>(j, "sa_sampling", s.sa_sampling);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseStoringFactors(j, s.storing_factor);
  return s;
}

inline PDLSweep ParsePDL(const nlohmann::json& j) {
  PDLSweep s;
  s.codec = ParseEnumList<axes::PDLVariant>(j, "codec", s.codec);
  s.get_doc = ParseEnumList<axes::GetDocEnum>(j, "get_doc", s.get_doc);
  s.policy = ParseEnumList<axes::PDLStoragePolicy>(j, "policy", s.policy);
  s.gcda_slp = ParseEnumList<axes::GCDASLPVariant>(j, "gcda_slp", s.gcda_slp);
  s.bare_slp = ParseEnumList<axes::BareSLPVariant>(j, "bare_slp", s.bare_slp);
  s.dgcda_slp = ParseEnumList<axes::DGCDASLPVariant>(j, "dgcda_slp", s.dgcda_slp);
  s.sa_sampling = ParseList<std::size_t>(j, "sa_sampling", s.sa_sampling);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseStoringFactors(j, s.storing_factor);
  s.backend_block_size = ParseList<std::uint32_t>(j, "backend_block_size", s.backend_block_size);
  return s;
}

inline GCDABackendSweep ParseGCDABackend(const nlohmann::json& j) {
  GCDABackendSweep s;
  s.get_doc = ParseEnumList<axes::GetDocEnum>(j, "get_doc", s.get_doc);
  s.sets = ParseEnumList<axes::GCDASetsCodec>(j, "sets", s.sets);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseStoringFactors(j, s.storing_factor);
  return s;
}

inline BruteSweep ParseBrute(const nlohmann::json& j) {
  BruteSweep s;
  s.kind = ParseListWith<BruteSweep::Kind>(j, "kind", s.kind, ParseBruteKind);
  s.sampling_size = ParseList<std::size_t>(j, "sampling_size", s.sampling_size);
  return s;
}

inline FamilySweep ParseFamilyBlock(const nlohmann::json& j, int t_version) {
  if (!j.contains("family") || !j.at("family").is_string()) {
    throw std::invalid_argument("sweep block: missing or non-string 'family'");
  }
  const auto family = j.at("family").get<std::string>();
  if (family == "gcda")    return ParseGCDA(j);
  if (family == "dgcda")   return ParseDGCDA(j);
  if (family == "slp_ns")  return ParseSLPNS(j);
  if (family == "rmq")     return ParseRMQ(j, t_version);
  if (family == "pdl")     return ParsePDL(j);
  if (family == "brute")   return ParseBrute(j);
  if (family == "gcda_backend") return ParseGCDABackend(j);
  throw std::invalid_argument("sweep block: unknown family '" + family + "'");
}

inline Dataset ParseDataset(const nlohmann::json& j) {
  Dataset d;
  if (j.contains("dir")) d.dir = j.at("dir").get<std::string>();
  if (j.contains("name")) d.name = j.at("name").get<std::string>();
  if (j.contains("patterns")) d.patterns = j.at("patterns").get<std::string>();
  if (j.contains("data_width")) d.data_width = j.at("data_width").get<std::uint8_t>();
  if (j.contains("sa_algo")) d.sa_algo = j.at("sa_algo").get<std::string>();
  if (j.contains("doc_delim")) d.doc_delim = j.at("doc_delim").get<int>();
  return d;
}

inline RePair ParseRePair(const nlohmann::json& j) {
  RePair r;
  if (j.contains("variant")) {
    r.variant = j.at("variant").get<std::string>();
    if (r.variant != "auto" && r.variant != "bal32" && r.variant != "bal64") {
      throw std::invalid_argument("spec.repair.variant: must be 'auto', 'bal32' or 'bal64', got '" +
                                  r.variant + "'");
    }
  }
  if (j.contains("mb")) r.mb = j.at("mb").get<std::uint64_t>();
  if (j.contains("mb_max")) r.mb_max = j.at("mb_max").get<std::uint64_t>();
  return r;
}

inline Workload ParseWorkload(const nlohmann::json& j) {
  Workload w;
  if (j.contains("reps")) w.reps = j.at("reps").get<int>();
  if (j.contains("min_time")) w.min_time = j.at("min_time").get<double>();
  if (j.contains("report_stats")) w.report_stats = j.at("report_stats").get<bool>();
  if (j.contains("print_result")) w.print_result = j.at("print_result").get<bool>();
  if (j.contains("rebuild")) w.rebuild = j.at("rebuild").get<bool>();
  if (j.contains("memory_trace")) w.memory_trace = j.at("memory_trace").get<bool>();
  if (j.contains("pattern_code")) w.pattern_code = j.at("pattern_code").get<std::string>();
  if (j.contains("pattern_delim")) w.pattern_delim = j.at("pattern_delim").get<int>();
  if (j.contains("task")) {
    const auto t = j.at("task").get<std::string>();
    if (t == "listing") w.task = Task::Listing;
    else if (t == "frequency") w.task = Task::Frequency;
    else throw std::invalid_argument("spec.workload.task: must be 'listing' or 'frequency', got '" + t + "'");
  }
  return w;
}

}  // namespace detail

inline Spec ParseSpec(const nlohmann::json& j) {
  Spec s;
  if (j.contains("spec_version")) {
    s.spec_version = j.at("spec_version").get<int>();
    if (s.spec_version < 1 || s.spec_version > detail::kSpecVersion)
      throw std::invalid_argument("spec_version: 1 (2026-09 names) or " + std::to_string(detail::kSpecVersion) +
                                  " (current), got " + std::to_string(s.spec_version));
  }
  if (j.contains("mode")) s.mode = detail::ParseMode(j.at("mode").get<std::string>());
  if (j.contains("dataset")) s.dataset = detail::ParseDataset(j.at("dataset"));
  if (j.contains("workload")) s.workload = detail::ParseWorkload(j.at("workload"));
  if (j.contains("repair")) s.repair = detail::ParseRePair(j.at("repair"));
  if (j.contains("sweep")) {
    for (const auto& block : j.at("sweep")) {
      s.sweep.push_back(detail::ParseFamilyBlock(block, s.spec_version));
    }
  }
  if (s.dataset.dir.empty() || s.dataset.name.empty()) {
    throw std::invalid_argument("spec.dataset: 'dir' and 'name' are required");
  }
  if (s.mode == Mode::Query && s.dataset.patterns.empty()) {
    throw std::invalid_argument("spec.dataset.patterns: required for query mode");
  }
  return s;
}

inline Spec ParseSpecFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open spec file: " + path);
  nlohmann::json j;
  in >> j;
  return ParseSpec(j);
}

}  // namespace bench::spec
