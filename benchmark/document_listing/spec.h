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
//     "mode":     "query" | "construct",
//     "dataset":  { "dir": "...", "name": "...", "patterns": "...",
//                   "data_width": 8, "sa_algo": "SDSL_LIBDIVSUFSORT",
//                   "doc_delim": 0 },
//     "workload": { "reps": 10, "min_time": 0.0,
//                   "report_stats": false, "print_result": false },
//     "sweep": [
//       { "family": "gcda",  "tslp": ["Light"], "block_size": [512],
//                            "storing_factor": [4] },
//       { "family": "dgcda", "tslp": ["Default","OTF"], "block_size": [512],
//                            "storing_factor": [4] },
//       { "family": "rmq",   "core": ["sada","ilcp"],
//                            "get_doc": ["da","slp"], "gcda_slp": ["Light"],
//                            "block_size": [512], "storing_factor": [4] },
//       { "family": "pdl",   "codec": ["Plain"], "get_doc": ["DA"],
//                            "policy": ["OccurrenceWeighted"],
//                            "block_size": [512], "storing_factor": [4] },
//       { "family": "slp_ns","tslp": ["IV"] },
//       { "family": "brute", "kind": ["r-index"], "sampling_size": [16] }
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
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <nlohmann/json.hpp>

#include "axes.h"
#include "enum_traits.h"
#include "factories/rmq.h"  // for rmq::CoreKind

namespace bench::spec {

enum class Mode { Query, Construct };

struct Dataset {
  std::string dir;
  std::string name;
  std::string patterns;  // Query mode only.
  std::uint8_t data_width = 8;
  std::string sa_algo = "SDSL_LIBDIVSUFSORT";
  int doc_delim = 0;
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
  // Default false (warm-cache fast-path, identical to bm_build_items.cpp).
  bool rebuild = false;
  // Construct mode only: write construction-<bench>.{html,json} memory-monitor
  // traces and sizes-<bench>.json size sidecars next to the working directory
  // for every benchmark cell — matches bm_build_items' always-on behaviour.
  bool memory_trace = false;
  // Pattern coding for query mode: "PLAIN" | "BASE64".
  std::string pattern_code = "PLAIN";
  int pattern_delim = '\n';
};

struct GCDASweep {
  std::vector<axes::GCDASLPVariant> tslp{axes::GCDASLPVariant::Light};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
};

struct DGCDASweep {
  std::vector<axes::DGCDASLPVariant> tslp{axes::DGCDASLPVariant::Default};
  std::vector<std::uint32_t> block_size{512};
  std::vector<float> storing_factor{4.0f};
};

struct SLPNSSweep {
  std::vector<axes::BareSLPVariant> tslp{axes::BareSLPVariant::IV};
};

struct RMQSweep {
  std::vector<factories::rmq::CoreKind> core{
      factories::rmq::CoreKind::SADA,
      factories::rmq::CoreKind::ILCP,
      factories::rmq::CoreKind::CILCP,
      factories::rmq::CoreKind::SADA_S,
      factories::rmq::CoreKind::ILCP_S,
      factories::rmq::CoreKind::CILCP_S};
  std::vector<axes::GetDocEnum> get_doc{axes::GetDocEnum::DA};
  std::vector<axes::GCDASLPVariant> gcda_slp{axes::GCDASLPVariant::Light};
  std::vector<axes::BareSLPVariant> bare_slp{axes::BareSLPVariant::IV};
  // DGCDA variant axis for get_doc=dslp (the differential get-doc backend).
  std::vector<axes::DGCDASLPVariant> dgcda_slp{axes::DGCDASLPVariant::Default};
  // TRunValues axis for the -S sub-family (ILCP-S / CILCP-S). Default DV
  // matches IlcpLikeSCore's default; other variants in {iv, dv, vv} fan
  // the -S cores out across the run-values container.
  std::vector<axes::RunValuesVariant> run_values{axes::RunValuesVariant::DV};
  // TPrevDoc axis for SADA-S. Default IV matches SadaSCore's default;
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
};

struct BruteSweep {
  enum class Kind { RIndex, SrIndex };
  std::vector<Kind> kind{Kind::RIndex};
  std::vector<std::size_t> sampling_size{8};  // Only used for sr-index.
};

using FamilySweep = std::variant<GCDASweep, DGCDASweep, SLPNSSweep,
                                  RMQSweep, PDLSweep, BruteSweep>;

struct Spec {
  Mode mode = Mode::Query;
  Dataset dataset;
  Workload workload;
  std::vector<FamilySweep> sweep;
};

//~~~~~~~  Parser  ~~~~~~~

namespace detail {

inline Mode ParseMode(const std::string& s) {
  if (s == "query")     return Mode::Query;
  if (s == "construct") return Mode::Construct;
  throw std::invalid_argument("spec.mode: must be 'query' or 'construct', got '" + s + "'");
}

inline factories::rmq::CoreKind ParseRmqCore(const std::string& s) {
  if (s == "sada")    return factories::rmq::CoreKind::SADA;
  if (s == "ilcp")    return factories::rmq::CoreKind::ILCP;
  if (s == "cilcp")   return factories::rmq::CoreKind::CILCP;
  if (s == "sada-s")  return factories::rmq::CoreKind::SADA_S;
  if (s == "ilcp-s")  return factories::rmq::CoreKind::ILCP_S;
  if (s == "cilcp-s") return factories::rmq::CoreKind::CILCP_S;
  throw std::invalid_argument("rmq.core: must be 'sada' / 'ilcp' / 'cilcp' / 'sada-s' / 'ilcp-s' / 'cilcp-s', got '" + s + "'");
}

inline BruteSweep::Kind ParseBruteKind(const std::string& s) {
  if (s == "r-index" || s == "rindex")    return BruteSweep::Kind::RIndex;
  if (s == "sr-index" || s == "srindex")  return BruteSweep::Kind::SrIndex;
  throw std::invalid_argument("brute.kind: must be 'r-index' or 'sr-index', got '" + s + "'");
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

inline GCDASweep ParseGCDA(const nlohmann::json& j) {
  GCDASweep s;
  s.tslp = ParseEnumList<axes::GCDASLPVariant>(j, "tslp", s.tslp);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseList<float>(j, "storing_factor", s.storing_factor);
  return s;
}

inline DGCDASweep ParseDGCDA(const nlohmann::json& j) {
  DGCDASweep s;
  s.tslp = ParseEnumList<axes::DGCDASLPVariant>(j, "tslp", s.tslp);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseList<float>(j, "storing_factor", s.storing_factor);
  return s;
}

inline SLPNSSweep ParseSLPNS(const nlohmann::json& j) {
  SLPNSSweep s;
  s.tslp = ParseEnumList<axes::BareSLPVariant>(j, "tslp", s.tslp);
  return s;
}

inline RMQSweep ParseRMQ(const nlohmann::json& j) {
  RMQSweep s;
  s.core = ParseListWith<factories::rmq::CoreKind>(j, "core", s.core, ParseRmqCore);
  s.get_doc = ParseEnumList<axes::GetDocEnum>(j, "get_doc", s.get_doc);
  s.gcda_slp = ParseEnumList<axes::GCDASLPVariant>(j, "gcda_slp", s.gcda_slp);
  s.bare_slp = ParseEnumList<axes::BareSLPVariant>(j, "bare_slp", s.bare_slp);
  s.dgcda_slp = ParseEnumList<axes::DGCDASLPVariant>(j, "dgcda_slp", s.dgcda_slp);
  s.run_values = ParseEnumList<axes::RunValuesVariant>(j, "run_values", s.run_values);
  s.prev_doc = ParseEnumList<axes::PrevDocVariant>(j, "prev_doc", s.prev_doc);
  s.sa_sampling = ParseList<std::size_t>(j, "sa_sampling", s.sa_sampling);
  s.block_size = ParseList<std::uint32_t>(j, "block_size", s.block_size);
  s.storing_factor = ParseList<float>(j, "storing_factor", s.storing_factor);
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
  s.storing_factor = ParseList<float>(j, "storing_factor", s.storing_factor);
  return s;
}

inline BruteSweep ParseBrute(const nlohmann::json& j) {
  BruteSweep s;
  s.kind = ParseListWith<BruteSweep::Kind>(j, "kind", s.kind, ParseBruteKind);
  s.sampling_size = ParseList<std::size_t>(j, "sampling_size", s.sampling_size);
  return s;
}

inline FamilySweep ParseFamilyBlock(const nlohmann::json& j) {
  if (!j.contains("family") || !j.at("family").is_string()) {
    throw std::invalid_argument("sweep block: missing or non-string 'family'");
  }
  const auto family = j.at("family").get<std::string>();
  if (family == "gcda")    return ParseGCDA(j);
  if (family == "dgcda")   return ParseDGCDA(j);
  if (family == "slp_ns")  return ParseSLPNS(j);
  if (family == "rmq")     return ParseRMQ(j);
  if (family == "pdl")     return ParsePDL(j);
  if (family == "brute")   return ParseBrute(j);
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
  return w;
}

}  // namespace detail

inline Spec ParseSpec(const nlohmann::json& j) {
  Spec s;
  if (j.contains("mode")) s.mode = detail::ParseMode(j.at("mode").get<std::string>());
  if (j.contains("dataset")) s.dataset = detail::ParseDataset(j.at("dataset"));
  if (j.contains("workload")) s.workload = detail::ParseWorkload(j.at("workload"));
  if (j.contains("sweep")) {
    for (const auto& block : j.at("sweep")) {
      s.sweep.push_back(detail::ParseFamilyBlock(block));
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
