//
// Generic parse / name helper for the doc-list benchmark axis enums.
//
// Each EnumTraits<E> specialisation knows how to convert a comma-separated
// CLI flag value into a std::vector<E> and how to render an E back to a
// human-readable string. Replaces the per-axis Parse*/*Name helper pairs
// that used to live duplicated in bm_query_doc_list.cpp and bm_build_items.cpp.
//

#pragma once

#include <algorithm>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "dret/pdl/storage_policy.h"

#include "axes.h"

namespace bench::axes {

template <typename E>
struct EnumTraits;  // primary undefined; specialisations below.

// Shared CSV-parse helper. Walks the comma-separated value, looks each token
// up via TryParse(token) -> std::optional<E>, throws on unknown tokens, falls
// back to a default if the result is empty.
template <typename E>
std::vector<E> ParseCSV(std::string_view t_value) {
  std::vector<E> out;
  std::stringstream ss{std::string(t_value)};
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item.empty()) continue;
    auto parsed = EnumTraits<E>::TryParse(item);
    if (!parsed) {
      throw std::invalid_argument(std::string(EnumTraits<E>::flag_name())
                                  + ": unknown value '" + item + "'");
    }
    out.push_back(*parsed);
  }
  if (out.empty() && EnumTraits<E>::has_default()) {
    out.push_back(EnumTraits<E>::default_value());
  }
  return out;
}

template <typename E>
bool Contains(const std::vector<E>& t_vec, E t_value) {
  return std::find(t_vec.begin(), t_vec.end(), t_value) != t_vec.end();
}

//~~~~~~~  Specialisations  ~~~~~~~

template <>
struct EnumTraits<GetDocEnum> {
  static constexpr const char* flag_name() { return "rmq_get_doc_variants"; }
  static constexpr bool has_default() { return true; }
  static constexpr GetDocEnum default_value() { return GetDocEnum::DA; }

  static std::optional<GetDocEnum> TryParse(std::string_view s) {
    if (s == "da")     return GetDocEnum::DA;
    if (s == "slp")    return GetDocEnum::SLP;
    if (s == "slp_ns") return GetDocEnum::SLP_NS;
    if (s == "dslp")   return GetDocEnum::DSLP;
    return std::nullopt;
  }

  static const char* Name(GetDocEnum v) {
    switch (v) {
      case GetDocEnum::DA:     return "DA";
      case GetDocEnum::SLP:    return "SLP";
      case GetDocEnum::SLP_NS: return "SLP-NS";
      case GetDocEnum::DSLP:   return "DSLP";
    }
    return "UNKNOWN";
  }
};

template <>
struct EnumTraits<GCDASLPVariant> {
  static constexpr const char* flag_name() { return "gcda_slp_variants"; }
  static constexpr bool has_default() { return true; }
  static constexpr GCDASLPVariant default_value() { return GCDASLPVariant::Default; }

  static std::optional<GCDASLPVariant> TryParse(std::string_view s) {
    if (s == "default")       return GCDASLPVariant::Default;
    if (s == "compact_bp")    return GCDASLPVariant::CompactBP;
    if (s == "compact_louds") return GCDASLPVariant::CompactLOUDS;
    if (s == "cslp")          return GCDASLPVariant::CSLP;
    return std::nullopt;
  }

  static const char* Name(GCDASLPVariant v) {
    switch (v) {
      case GCDASLPVariant::Default:      return "Default";
      case GCDASLPVariant::CompactBP:    return "CompactBP";
      case GCDASLPVariant::CompactLOUDS: return "CompactLOUDS";
      case GCDASLPVariant::CSLP:         return "CSLP";
    }
    return "UNKNOWN";
  }
};

template <>
struct EnumTraits<BareSLPVariant> {
  static constexpr const char* flag_name() { return "bare_slp_variants"; }
  static constexpr bool has_default() { return true; }
  static constexpr BareSLPVariant default_value() { return BareSLPVariant::Default; }

  static std::optional<BareSLPVariant> TryParse(std::string_view s) {
    if (s == "default") return BareSLPVariant::Default;
    if (s == "raw")     return BareSLPVariant::Raw;
    if (s == "dv")      return BareSLPVariant::DV;
    if (s == "vv")      return BareSLPVariant::VV;
    return std::nullopt;
  }

  static const char* Name(BareSLPVariant v) {
    switch (v) {
      case BareSLPVariant::Default: return "Default";
      case BareSLPVariant::Raw:     return "Raw";
      case BareSLPVariant::DV:      return "DV";
      case BareSLPVariant::VV:      return "VV";
    }
    return "UNKNOWN";
  }
};

template <>
struct EnumTraits<PDLVariant> {
  static constexpr const char* flag_name() { return "pdl_variants"; }
  // No default: when --pdl_variants is empty, PDL is disabled entirely.
  static constexpr bool has_default() { return false; }
  static constexpr PDLVariant default_value() { return PDLVariant::Plain; }  // unused

  static std::optional<PDLVariant> TryParse(std::string_view s) {
    if (s == "plain") return PDLVariant::Plain;
    if (s == "rp")    return PDLVariant::RP;
    if (s == "bc")    return PDLVariant::BC;
    return std::nullopt;
  }

  static const char* Name(PDLVariant v) {
    switch (v) {
      case PDLVariant::Plain: return "Plain";
      case PDLVariant::RP:    return "RP";
      case PDLVariant::BC:    return "BC";
    }
    return "UNKNOWN";
  }
};

template <>
struct EnumTraits<PDLStoragePolicy> {
  static constexpr const char* flag_name() { return "pdl_storage_policy"; }
  static constexpr bool has_default() { return true; }
  static constexpr PDLStoragePolicy default_value() { return PDLStoragePolicy::OccurrenceWeighted; }

  static std::optional<PDLStoragePolicy> TryParse(std::string_view s) {
    if (s == "occurrence_weighted") return PDLStoragePolicy::OccurrenceWeighted;
    if (s == "all_internal")        return PDLStoragePolicy::StoreAllInternal;
    if (s == "leaves_only")         return PDLStoragePolicy::LeavesOnly;
    return std::nullopt;
  }

  static const char* Name(PDLStoragePolicy v) {
    switch (v) {
      case PDLStoragePolicy::OccurrenceWeighted: return "OccurrenceWeighted";
      case PDLStoragePolicy::StoreAllInternal:   return "StoreAllInternal";
      case PDLStoragePolicy::LeavesOnly:         return "LeavesOnly";
    }
    return "UNKNOWN";
  }
};

// Map the factory-side PDLStoragePolicy enum to the library-side
// dret::pdl::StoragePolicy. Used in both binaries to pass the runtime
// argument into the PDL index constructor.
inline dret::pdl::StoragePolicy toPDLStoragePolicy(PDLStoragePolicy p) {
  switch (p) {
    case PDLStoragePolicy::StoreAllInternal:   return dret::pdl::StoragePolicy::StoreAllInternal;
    case PDLStoragePolicy::LeavesOnly:         return dret::pdl::StoragePolicy::LeavesOnly;
    case PDLStoragePolicy::OccurrenceWeighted:
    default:                                   return dret::pdl::StoragePolicy::OccurrenceWeighted;
  }
}

}  // namespace bench::axes
