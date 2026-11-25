//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Benchmark names: FAMILY[k=v|k=v], the keys in the order pushed, an empty value left out. Callers
// follow the canonical order codec · get-doc · slp · span-length · slp-container · policy · tree ·
// block-size · storing-factor · prev-doc · run-values · sampling. One name per cell, the same in query
// and construct mode, and the analysis scripts read them (experiments/analyze.py).
//

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "dret/config.h"

#include "../axes.h"

namespace bench::names {

using KV = std::vector<std::pair<std::string, std::string>>;

inline std::string KeyedName(const std::string& t_family, const KV& t_kvs) {
  std::string s = t_family;
  bool first = true;
  for (const auto& [k, v] : t_kvs) {
    if (v.empty()) continue;
    s += first ? '[' : '|';
    s += k;
    s += '=';
    s += v;
    first = false;
  }
  if (!first) s += ']';
  return s;
}

inline std::string IntStr(std::uint64_t t_x) {
  return std::to_string(t_x);
}

inline std::string SfStr(float t_sf) {
  return dret::IsInfiniteStoringFactor(t_sf) ? std::string("inf") : std::to_string(static_cast<int>(t_sf));
}

inline const char* CoreName(axes::CoreKind t_core) {
  switch (t_core) {
    case axes::CoreKind::SADA:
      return "SADA";
    case axes::CoreKind::ILCP:
      return "ILCP";
    case axes::CoreKind::CILCP:
      return "CILCP";
    case axes::CoreKind::SADA_V:
      return "SADA-V";
    case axes::CoreKind::ILCP_V:
      return "ILCP-V";
    case axes::CoreKind::CILCP_V:
      return "CILCP-V";
  }
  return "UNKNOWN";
}

inline const char* SlpPlainValue(axes::GCDASLPVariant t_v) {
  switch (t_v) {
    case axes::GCDASLPVariant::CompactBP:
      return "compact-bp";
    case axes::GCDASLPVariant::CompactLOUDS:
      return "compact-louds";
    case axes::GCDASLPVariant::Combined:
      return "combined";
    case axes::GCDASLPVariant::Light:
    default:
      return "light";
  }
}

// A differential SLP variant decomposes into (span-length, slp-container).
inline std::pair<std::string, std::string> DiffParts(axes::DGCDASLPVariant t_v) {
  switch (t_v) {
    case axes::DGCDASLPVariant::OTF:
      return {"on-the-fly", "iv"};
    case axes::DGCDASLPVariant::CRL:
      return {"cached-root", "iv"};
    case axes::DGCDASLPVariant::EV:
      return {"full", "ev"};
    case axes::DGCDASLPVariant::DV:
      return {"full", "dv"};
    case axes::DGCDASLPVariant::VV:
      return {"full", "vv"};
    case axes::DGCDASLPVariant::Default:
    default:
      return {"full", "iv"};
  }
}

inline const char* BareValue(axes::BareSLPVariant t_v) {
  switch (t_v) {
    case axes::BareSLPVariant::Raw:
      return "raw";
    case axes::BareSLPVariant::DV:
      return "dv";
    case axes::BareSLPVariant::VV:
      return "vv";
    case axes::BareSLPVariant::IV:
    default:
      return "iv";
  }
}

// The four bare-diff variants (base differential SLP, non-sampled).
inline bool IsBareDiff(axes::BareSLPVariant t_v) {
  return t_v == axes::BareSLPVariant::Diff || t_v == axes::BareSLPVariant::DiffEV || t_v == axes::BareSLPVariant::DiffDV
         || t_v == axes::BareSLPVariant::DiffVV;
}

// Differential int-container for a bare-diff variant (roots/span_sums/samples).
inline const char* DiffContainerValue(axes::BareSLPVariant t_v) {
  switch (t_v) {
    case axes::BareSLPVariant::DiffEV:
      return "ev";
    case axes::BareSLPVariant::DiffDV:
      return "dv";
    case axes::BareSLPVariant::DiffVV:
      return "vv";
    case axes::BareSLPVariant::Diff:
    default:
      return "iv";
  }
}

// The slp keys of a bare grammar: a differential one says slp=diff (BareValue maps every diff variant
// to "iv", so without it a differential grammar would take the plain iv grammar's name).
inline KV BareKeys(axes::BareSLPVariant t_v) {
  if (IsBareDiff(t_v)) return {{"slp", "diff"}, {"slp-container", DiffContainerValue(t_v)}};
  return {{"slp-container", BareValue(t_v)}};
}

// GCDA document-set codec in keyed form. Re-Pair, the family default, is the empty value, which
// KeyedName omits.
inline const char* SetsValue(axes::GCDASetsCodec t_v) {
  return t_v == axes::GCDASetsCodec::Plain ? "plain" : "";
}

inline const char* RunValuesValue(axes::RunValuesVariant t_v) {
  switch (t_v) {
    case axes::RunValuesVariant::DV:
      return "dv";
    case axes::RunValuesVariant::VV:
      return "vv";
    case axes::RunValuesVariant::IV:
    default:
      return "iv";
  }
}

inline const char* PrevDocValue(axes::PrevDocVariant t_v) {
  switch (t_v) {
    case axes::PrevDocVariant::DV:
      return "dv";
    case axes::PrevDocVariant::VV:
      return "vv";
    case axes::PrevDocVariant::IV:
    default:
      return "iv";
  }
}

inline const char* CodecValue(axes::PDLVariant t_v) {
  switch (t_v) {
    case axes::PDLVariant::RP:
      return "rp";
    case axes::PDLVariant::BC:
      return "bc";
    case axes::PDLVariant::Plain:
    default:
      return "plain";
  }
}

inline const char* PolicyValue(axes::PDLStoragePolicy t_v) {
  switch (t_v) {
    case axes::PDLStoragePolicy::LeavesOnly:
      return "leaves-only";
    case axes::PDLStoragePolicy::ListWeighted:
      return "list-weighted";
    case axes::PDLStoragePolicy::StoreAllInternal:
      return "store-all";
    case axes::PDLStoragePolicy::OccurrenceWeighted:
    default:
      return "occ-weighted";
  }
}

// get-doc backend kind in keyed form.
inline const char* GetDocValue(axes::GetDocEnum t_g) {
  switch (t_g) {
    case axes::GetDocEnum::SLP:
      return "gcda";
    case axes::GetDocEnum::SLP_NS:
      return "gcda-bare";
    case axes::GetDocEnum::DSLP:
      return "gcda-diff";
    case axes::GetDocEnum::SAPhiR:
      return "sa-phi-r";
    case axes::GetDocEnum::SAPhiSR:
      return "sa-phi-sr";
    case axes::GetDocEnum::RLCSA:
      return "rlcsa";
    case axes::GetDocEnum::DA:
    default:
      return "da";
  }
}

}  // namespace bench::names
