//
// Axis enums of the doc-list benchmark (bm_doc_list).
//
// These enums identify the variation axes across the doc-list index family
// (which family / which SLP / which GetDoc policy / which PDL codec / which
// PDL storage policy), shared by the spec parser, the factories and the
// driver. Factory<>::GetDocEnum etc. are re-exported as type aliases.
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
  BRUTE_DA,     // dret::DocListIdxBruteBackend over the plain DA (counts with the r-index)
  BRUTE_RLCSA,  // dret::DocListIdxBruteBackend over the RLCSA (counts with the RLCSA)
  COUNT_RI,     // the r-index count alone (factories::brute::IdxCount): no document is reported
  GCDA,
  DGCDA,    // variant axis selects Default / OTF / CRL / EV / DV / VV (see DGCDASLPVariant)
  // These names, and the benchmark labels built from them, are the 2026-11
  // vocabulary: the plain names are the published algorithms (marking, no
  // stored values) and the -V names dret's value-rule variants. Producers stamp
  // context.core_naming = "2026-11". Under the 2026-09 names (results before
  // then) the marking cores were SADA-L / ILCP-L / CILCP-L and the value cores
  // had the plain names.
  SADA,     // -> SadaCore:  RMinQ on prev_doc, marker-based stop
  ILCP,     // -> IlcpCore:  RMinQ on backward-ILCP runs, marker-based stop
  CILCP,    // -> CilcpCore: CMR20's document-aware merged runs, marker-based stop
  SLP_NS,   // dret::DocListIdxSLP (Brute-GCDA) — brute force over a GCDA grammar
  PDL,      // Precomputed Document Listing — variant axis selects Plain/RP/BC
  // dret's value-rule variants: same runs as the cores above, plus the stored
  // array their value-based recursion stop consults.
  SADA_V,   // -> SadaValueCore
  ILCP_V,   // -> IlcpValueCore
  CILCP_V,  // -> CilcpValueCore: the SAME CMR20 partition as CILCP,
            // built by one shared routine; the two differ only in whether the
            // per-run minima are stored, not in where the runs fall.
  GCDA_BACKEND,  // GCDA's sampled tree and node lists over a DA backend (DA, SA-Phi, RLCSA)
};

// PDL stored-set codec axis. The class template differs per value
// (DocListIdxPDLPlain / DocListIdxPDLRP / DocListIdxPDLBC).
enum class PDLVariant {
  Plain,
  RP,
  BC,
};

// How GCDA stores the document sets of its sampled nodes: Re-Pair compressed
// (grammar::GCChunks, the family default) or plain -- the sorted lists
// themselves, bit-packed (grammar::Chunks of int_vector), as PDL's plain codec
// stores them. Applies to every TSLP of the GCDA and DGCDA families.
enum class GCDASetsCodec {
  RP,
  Plain,
};

// PDL tree-construction storage policy. Threaded into the constructor as a
// runtime argument; baked into the cache-key prefix.
enum class PDLStoragePolicy {
  OccurrenceWeighted,
  StoreAllInternal,
  LeavesOnly,
  ListWeighted,
};

// RMQ raw-range get-doc policy. The values DA / SLP / DSLP apply to both
// RMQ-listing cores (SADA / ILCP / CILCP) and the PDL family. SLP_NS is only
// meaningful for the RMQ family (it consumes the bare grammar::SLP<> cache).
//
// SAPhiR / SAPhiSR are RLCSA-style backings (PDL paper baseline): compute
// DA[i] = rank1(doc_ends, SA[i]) where SA[i] is recovered via Phi-walks from
// run-end samples (r-index, dense) or subsampled run-end samples (sr-index,
// with subsample_rate as the runtime sa_sampling sub-axis). Shared class
// (rmq::GetDocSAPhi) consumed by both RMQ cores and PDL via the same
// dispatch pattern as DA/SLP/SLP_NS/DSLP. See docs/pdl_rlcsa_baseline_plan.md.
enum class GetDocEnum {
  DA,
  SLP,
  SLP_NS,
  DSLP,
  SAPhiR,   // r-index-backed (no sampling param; dense, RLCSA-analogue)
  SAPhiSR,  // sr-index-backed (carries sa_sampling sub-axis; subsampled)
  RLCSA,    // paper-faithful: batched CSA::RLCSA::locate(range) + getSequenceForPosition
};

// GCDA's TSLP choice. The RMQ-SLP path also consults this so SADA / ILCP /
// CILCP-SLP reuse the SLP cache file built for the matching GCDA variant.
// `Light` is the canonical GCDA TSLP (grammar::LightSLP<...>) — the family
// default; the other three are alternative grammar representations.
enum class GCDASLPVariant {
  Light,         // grammar::LightSLP<...> — canonical GCDA TSLP; the family default
  CompactBP,     // grammar::CompactBPSLP<>
  CompactLOUDS,  // grammar::CompactLOUDSSLP<>
  Combined,          // grammar::CombinedSLPWithUnitCover<> (sampled-ondemand)
};

// Bare-SLP container choice for the SLP-NS family. enc_vector<> is excluded:
// rule pairs and span lengths are non-monotonic. The four values name the
// inner container in (IV=int_vector / Raw=library default vector<uint32_t> /
// DV=dac_vector / VV=vlc_vector).
enum class BareSLPVariant {
  IV,       // grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>> — the family default
  Raw,      // grammar::SLP<> — library defaults (std::vector<uint32_t>); DRL-equivalent
  DV,       // grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>
  VV,       // grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>
  // bare-diff: base dret::DifferentialSLP<> (non-sampled differential SLP, no
  // GCChunks), with the differential int-container axis on roots/span_sums/
  // samples — Diff=int_vector(iv), DiffEV=enc_vector, DiffDV=dac_vector,
  // DiffVV=vlc_vector (mirrors the DGCDA sampled-diff EV/DV/VV variants). Only
  // the SLP-NS *listing* family consumes these; RMQ/PDL bare get-doc backends
  // are plain-only and ignore them.
  Diff,
  DiffEV,
  DiffDV,
  DiffVV,
};

// TRunValues container choice for the value-rule doc-listing cores
// (IndexEnum::ILCP_V / CILCP). The persisted per-run min(VILCP) array
// can be encoded fixed-width and bit-compressed (IV) or with variable-
// length per-element codes (DV / VV).
enum class RunValuesVariant {
  IV,       // sdsl::int_vector<>  — fixed-width, sdsl::util::bit_compress
  DV,       // sdsl::dac_vector<>  — direct access codes (family default)
  VV,       // sdsl::vlc_vector<>  — variable-length codes
};

// TPrevDoc container choice for the Sadakane-style SADA family
// (IndexEnum::SADA_V). The
// persisted prev_doc array stores per-SA-position previous-occurrence
// SA positions — values cover [0, n) roughly uniformly, so the natural
// default is fixed-width int_vector (DAC / VLC have little to compress
// in a uniform distribution); the axis exists so users can verify that
// empirically without recompiling.
enum class PrevDocVariant {
  IV,       // sdsl::int_vector<>  — fixed-width, sdsl::util::bit_compress (family default)
  DV,       // sdsl::dac_vector<>
  VV,       // sdsl::vlc_vector<>
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
