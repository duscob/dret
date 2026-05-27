//
// SA-Phi document-array lookup policy — the RLCSA-style backing for the PDL
// paper baseline (Cobas/Mäkinen/Rossi SPIRE 2020). Computes
//
//     DA[i] = rank1(doc_ends, SA[i])
//
// where SA[i] is recovered via Phi-walks from r-index run-end samples (or
// sr-index subsampled run-end samples when TLocateIdx is an SrIndex
// specialisation — deferred; see docs/pdl_rlcsa_baseline_plan.md).
//
// Same interface as the rest of the rmq::GetDoc* family
// (rmq_get_doc_policies.h) so both PDL (via PDLRawRangePolicy) and the RMQ
// cores can consume it through their existing TGetDoc template parameter.
//
// The class is templated on TLocateIdx so the same implementation serves
// both backings:
//   - sri::RIndex            — dense run-end samples (the RLCSA-analogue).
//   - sri::SrIndexValidArea  — subsampled (subsample_rate runtime knob);
//                              currently treated as RIndex; see plan §B.2.
//
// SA[i] algorithm:
//   1. Find the BWT run r containing position i  →  bwt_rle.run_of_position(i).
//   2. Look up the run-end position             →  [_, run_end] = run_range(r).
//   3. Read the run-end SA sample                →  samples[r] + 1  =  SA[run_end].
//      (The +1 mirrors the convention used by sri::RIndex's toehold; see
//      r_index.h:304.)
//   4. Walk Phi backward (run_end - i) times:
//        sa = SA[run_end];  for k in (run_end..=i+1): sa = phi(sa).first;
//      yields SA[i].
//
// All loaded items are shared (via the typed-cache + std::any storage map)
// with `BruteRI` / `BruteSRI` / `GetDocBv` — no extra files on disk.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>

#include "sr-index/construct.h"           // sri::construct for r-index
#include "sr-index/construct_base.h"      // sri::conf::KEY_BWT_*
#include "sr-index/phi.h"                 // sri::buildPhiBackward
#include "sr-index/r_index.h"             // sri::RIndex (default TLocateIdx)
#include "sr-index/rle_string.hpp"        // sri::RLEStringS
#include "sr-index/sampling.h"            // sri::SampleValidatorDefault
#include "sr-index/sequence_ops.h"        // sri::CircularPredecessor, RandomAccessFor*

#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret::rmq {

// Document-array lookup over SA[i] + doc-ends bitvector. Templated on
// TLocateIdx so the same class serves both r-index and sr-index backings —
// matching the BruteRI / BruteSRI split, reused here as a TGetDoc policy for
// both RMQ cores and PDL.
template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TLocateIdx = sri::RIndex<TStorage, dret::Alphabet<t_width>>,
          typename TBvDocEnds = sdsl::sd_vector<>>
class GetDocSAPhi : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  // For TLocateIdx = SrIndex variants, t_sa_sampling is the subsample rate;
  // for TLocateIdx = RIndex it is ignored (r-index has no sampling knob).
  // SR support is currently deferred; treat as RIndex.
  explicit GetDocSAPhi(const TStorage& t_storage, std::size_t t_sa_sampling = 0)
      : Base(t_storage), sa_sampling_(t_sa_sampling) {}

  GetDocSAPhi() = default;

  std::size_t sa_sampling() const { return sa_sampling_; }

  // Single-position lookup. Returns DA[i] = rank1(doc_ends, SA[i]).
  std::size_t operator()(std::size_t t_i) const {
    return doc_at_(sa_at_(t_i));
  }

  // Range expansion [b, e). Mirrors the existing rmq::GetDoc* API so PDL's
  // PDLRawRangePolicy can wrap this verbatim.
  template <typename TReport>
  void operator()(std::size_t t_b, std::size_t t_e, TReport& t_report) const {
    for (std::size_t i = t_b; i < t_e; ++i) {
      t_report((*this)(i));
    }
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;
    // The locate index + doc_ends bitvector live in shared storage with
    // BruteRI / BruteSRI / GetDocBv (typed-cache hash sharing). Per-instance
    // serialised state is just the sampling-rate scalar.
    written_bytes += sdsl::write_member(sa_sampling_, out, child, "sa_sampling");
    return written_bytes;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    // Shared caches (bwt_rle, samples, marks, mark_to_sample, doc_ends) are
    // owned upstream by the brute baselines + GetDocBv — accounted there.
    // GetDocSAPhi's marginal cost is just the scalar sa_sampling_.
    return r;
  }

 protected:
  using typename Base::TSource;
  using TBwtRLE = sri::RLEStringS<t_width>;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace sri::conf;

    // r-index data structures (shared with BruteRI's loaded instance; loadItem
    // returns the typed-cache ref_wrapper). These are NOT type-hashed on the
    // disk cache key (they're the standard sr-index sidecar files), so we
    // call loadItem with add_type_hash=false (the default).
    auto cref_bwt        = this->template loadItem<TBwtRLE>(KEY_BWT_RLE, t_source);
    auto cref_samples    = this->template loadItem<sdsl::int_vector<>>(KEY_BWT_RUN_LAST_TEXT_POS, t_source);
    auto cref_m2s        = this->template loadItem<sdsl::int_vector<>>(KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_TO_LAST_IDX, t_source);

    // Marks bitvector + rank/select supports (use loadItemPtr-backed BVRank/Select
    // to guarantee stable pointers across std::any reallocations).
    this->template loadItem<sdsl::sd_vector<>>(KEY_BWT_RUN_FIRST_TEXT_POS, t_source, true);
    auto bv_marks_rank   = this->template loadBVRank<sdsl::sd_vector<>>(KEY_BWT_RUN_FIRST_TEXT_POS, t_source, true);
    auto bv_marks_select = this->template loadBVSelect<sdsl::sd_vector<>>(KEY_BWT_RUN_FIRST_TEXT_POS, t_source, true);

    const std::size_t n = cref_bwt.get().size();

    // Build the Phi-backward callable using the same pieces RIndex builds
    // internally (r_index.h:316-332).
    auto predecessor      = sri::CircularPredecessor(bv_marks_rank, bv_marks_select, n);
    auto get_sample_idx   = sri::RandomAccessForTwoContainersDefault(cref_m2s, /*default_validity=*/true);
    auto get_sample       = sri::RandomAccessForCRefContainer(cref_samples);
    auto sample_validator = sri::SampleValidatorDefault();
    auto phi              = sri::buildPhiBackward(predecessor, get_sample_idx, get_sample, sample_validator, n);

    // SA[i] = walk Phi from the run-end sample down to position i.
    // The public RLEString::rank(i, report) reports
    //   (rnk_in_bwt, c, run_rnk, run_start, run_end, symbol_run_rnk)
    // where run_end is exclusive (the position after the last position of
    // the run). The run-end sample at index run_rnk corresponds to
    // position run_end - 1.
    sa_at_ = [cref_bwt, cref_samples, phi, n](std::size_t i) -> std::size_t {
      std::size_t run_rnk = 0;
      std::size_t run_end_excl = 0;
      cref_bwt.get().rank(i, [&run_rnk, &run_end_excl](
          auto /*rnk*/, auto /*c*/, auto rr, auto /*rs*/, auto re, auto /*srr*/) {
        run_rnk = rr;
        run_end_excl = re;
      });
      std::size_t sa = (cref_samples.get()[run_rnk] + 1) % n;  // SA[run_end_excl - 1]
      std::size_t pos = run_end_excl - 1;
      while (pos > i) {
        sa = phi(sa).first;
        --pos;
      }
      return sa;
    };

    // doc_at_(sa) = rank1(doc_ends, sa). Same construction as GetDocBv
    // (doc_list_brute.h:124-132); shares the doc_ends cache via typed-cache.
    // The logical key kDocEnds maps to the actual cache file key
    // ("doc_end") via the JSON keys mapping.
    const auto key_doc_ends = t_keys[dret::conf::kDocEnds].template get<std::string>();
    this->template loadItem<TBvDocEnds>(key_doc_ends, t_source, true);
    auto doc_ends_rank = this->template loadBVRank<TBvDocEnds>(key_doc_ends, t_source, true);
    doc_at_ = [doc_ends_rank](std::size_t pos) -> std::size_t {
      return doc_ends_rank.get()(pos);
    };
  }

  std::function<std::size_t(std::size_t)> sa_at_;   // i   -> SA[i]
  std::function<std::size_t(std::size_t)> doc_at_;  // pos -> rank1(doc_ends, pos)
  std::size_t sa_sampling_ = 0;
};

// Free-function construct() — same shape as the other rmq::GetDoc* construct
// overloads in rmq_get_doc_policies.h. Builds (or short-circuits when warm)
// the underlying r-index cache + the doc-ends bitvector.
template <typename TStorage, uint8_t t_width, typename TLocateIdx, typename TBvDocEnds>
void construct(GetDocSAPhi<TStorage, t_width, TLocateIdx, TBvDocEnds>& t_get_doc, Config& t_config) {
  using namespace dret::conf;

  // Text + doc-ends caches (shared with GetDocBv / brute / etc.). Built once
  // per collection; cheap no-op when already present.
  if (!cache_file_exists(t_config.keys[kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<t_width>(t_config);
  }
  if (!sdsl::cache_file_exists<TBvDocEnds>(t_config.keys[kDocEnds].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("DocEnds");
    ConstructDocEnd<t_width, TBvDocEnds>(t_config);
  }

  // r-index cache (shared with BruteRI). The sri::construct overload builds
  // the bwt_rle, samples (KEY_BWT_RUN_LAST_TEXT_POS), marks
  // (KEY_BWT_RUN_FIRST_TEXT_POS), and mark_to_sample
  // (KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_TO_LAST_IDX) sidecar files. Idempotent
  // when the caches already exist — cheap when warm.
  TLocateIdx tmp_idx(t_get_doc.storage());
  sri::construct(tmp_idx, t_config.data_path.string(), t_config);
}

}  // namespace dret::rmq
