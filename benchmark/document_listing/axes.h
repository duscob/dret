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
enum class IndexEnum {
  BRUTE_R_INDEX,
  BRUTE_SR_INDEX,
  GCDA,
  DGCDA,
  DGCDA_OTF,  // BasicSLPOnTheFlySpanLength — no stored lengths
  DGCDA_CRL,  // BasicSLPCachedRootSpanLengths — lengths cached for roots only
  DGCDA_EV,   // default SLP; TRoots/TSpanSums/TSamples = sdsl::enc_vector<>
  DGCDA_DV,   // default SLP; TRoots/TSpanSums/TSamples = sdsl::dac_vector<>
  DGCDA_VV,   // default SLP; TRoots/TSpanSums/TSamples = sdsl::vlc_vector<>
  SADA,       // RMinQ on prev_doc
  ILCP,       // RMinQ on backward-ILCP runs
  CILCP,      // RMinQ on doc-aware compressed backward-ILCP runs
  SLP_NS,     // Phase C: dret::DocListIdxSLP — non-sampled grammar::SLP<>
  PDL,        // Precomputed Document Listing — variant axis selects Plain/RP/BC
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

}  // namespace bench::axes
