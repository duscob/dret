//
// Shared axis enums for the doc-list benchmark binaries.
//
// These enums identify the variation axes across the doc-list index family
// (which family / which SLP / which GetDoc policy / which PDL codec / which
// PDL storage policy). They live at namespace scope so both bm_query_doc_list
// and bm_build_items can use the same types instead of redeclaring locally.
// Factory<>::GetDocEnum etc. are re-exported as type aliases for source
// compatibility with existing call sites.
//

#pragma once

namespace bench::axes {

// Which doc-list family / variant the benchmark is targeting.
//
// DGCDA's TSLP sub-axis (Default / OTF / CRL / EV / DV / VV) lives in
// DGCDASLPVariant below; consumers fan out across that axis just like GCDA fans
// across GCDASLPVariant.
enum class IndexEnum {
  BRUTE_R_INDEX,
  BRUTE_SR_INDEX,
  GCDA,
  DGCDA,    // variant axis selects Default / OTF / CRL / EV / DV / VV (see DGCDASLPVariant)
  SADA,     // RMinQ on prev_doc
  ILCP,     // RMinQ on backward-ILCP runs
  CILCP,    // RMinQ on doc-aware compressed backward-ILCP runs
  SLP_NS,   // Phase C: dret::DocListIdxSLP — non-sampled grammar::SLP<>
  PDL,      // Precomputed Document Listing — variant axis selects Plain/RP/BC
};

// PDL stored-set codec axis. The class template differs per value
// (DocListIdxPDLPlain / DocListIdxPDLRP / DocListIdxPDLBC).
enum class PDLVariant {
  Plain,
  RP,
  BC,
};

// PDL tree-construction storage policy. Threaded into the constructor as a
// runtime argument; baked into the cache-key prefix.
enum class PDLStoragePolicy {
  OccurrenceWeighted,
  StoreAllInternal,
  LeavesOnly,
};

// RMQ raw-range get-doc policy. The values DA / SLP / DSLP apply to both
// RMQ-listing cores (SADA / ILCP / CILCP) and the PDL family. SLP_NS is only
// meaningful for the RMQ family (it consumes the bare grammar::SLP<> cache).
enum class GetDocEnum {
  DA,
  SLP,
  SLP_NS,
  DSLP,
};

// GCDA's TSLP choice. The RMQ-SLP path also consults this so SADA / ILCP /
// CILCP-SLP reuse the SLP cache file built for the matching GCDA variant.
enum class GCDASLPVariant {
  Default,       // grammar::LightSLP<...> — the existing GCDA default
  CompactBP,     // grammar::CompactBPSLP<>
  CompactLOUDS,  // grammar::CompactLOUDSSLP<>
  CSLP,          // grammar::CombinedSLPWithUnitCover<> — Phase B
};

// Bare-SLP container choice for the SLP-NS family. enc_vector<> is excluded:
// rule pairs and span lengths are non-monotonic.
enum class BareSLPVariant {
  Default,  // grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>
  Raw,      // grammar::SLP<> — library defaults (std::vector<uint32_t>); DRL-equivalent
  DV,       // grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>
  VV,       // grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>
};

// DGCDA's TSLP choice. Default = the standard DifferentialLightSLP<>; OTF /
// CRL vary the BasicSLP span-length strategy; EV / DV / VV vary the inner
// int-vector container for roots / span_sums / samples.
enum class DGCDASLPVariant {
  Default,  // DifferentialLightSLP<> — the existing DGCDA default
  OTF,      // BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>
  CRL,      // BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>
  EV,       // grammar::SLP<>; TRoots/TSpanSums/TSamples = sdsl::enc_vector<>
  DV,       // grammar::SLP<>; TRoots/TSpanSums/TSamples = sdsl::dac_vector<>
  VV,       // grammar::SLP<>; TRoots/TSpanSums/TSamples = sdsl::vlc_vector<>
};

}  // namespace bench::axes
