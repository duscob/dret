//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/22/26.
//

#pragma once

#include <sdsl/bit_vectors.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <grammar/differential_slp.h>
#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>
#include <grammar/utility.h>

#include "dret/slp/basic_slp_span_length.h"
#include "dret/config.h"
#include "dret/size_report.h"

namespace dret {

template <typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
          typename TRoots = sdsl::int_vector<>,
          typename TSpanSums = sdsl::int_vector<>,
          typename TSamples = sdsl::int_vector<>,
          typename TSampleRootsPos = sdsl::enc_vector<>,
          typename TBV = sdsl::sd_vector<>>
class DifferentialSLP : public TSLP {
 public:
  using size_type = std::size_t;

  DifferentialSLP() = default;

  template <typename TOtherSLP,
            typename TOtherRoots,
            typename TOtherSpanSums,
            typename TOtherSamples,
            typename TOtherSampleRootsPos,
            typename TOtherBV,
            typename ActionSLPRules = grammar::NoAction,
            typename ActionSLPLengths = grammar::NoAction,
            typename ActionIntContainers = grammar::NoAction>
  DifferentialSLP(
      const DifferentialSLP<TOtherSLP, TOtherRoots, TOtherSpanSums, TOtherSamples, TOtherSampleRootsPos, TOtherBV>&
          other,
      ActionSLPRules&& action_slp_rules = grammar::NoAction(),
      ActionSLPLengths&& action_slp_lengths = grammar::NoAction(),
      ActionIntContainers&& action_int_containers = grammar::NoAction())
      : TSLP(static_cast<const TOtherSLP&>(other),
             std::forward<ActionSLPRules>(action_slp_rules),
             std::forward<ActionSLPLengths>(action_slp_lengths)),
        seq_size_{other.SeqSize()},
        diff_base_seq_{other.DiffBaseSeq()},
        diff_base_sums_{other.DiffBaseSums()} {
    grammar::Construct(roots_, other.GetRoots());
    action_int_containers(roots_);
    grammar::Construct(span_sums_, other.GetSpanSums());
    action_int_containers(span_sums_);
    grammar::Construct(samples_, other.GetSamples());
    action_int_containers(samples_);
    grammar::Construct(sample_roots_pos_, other.GetSampleRootsPos());
    action_int_containers(sample_roots_pos_);

    grammar::Construct(samples_pos_, other.GetSamplesPos());
    samples_pos_rank_ = typename TBV::rank_1_type(&samples_pos_);
    samples_pos_select_ = typename TBV::select_1_type(&samples_pos_);
  }

  const auto& GetRoots() const {
    return roots_;
  }

  const auto& GetSpanSums() const {
    return span_sums_;
  }

  const auto& GetSamples() const {
    return samples_;
  }

  const auto& GetSampleRootsPos() const {
    return sample_roots_pos_;
  }

  const auto& GetSamplesPos() const {
    return samples_pos_;
  }

  auto SeqSize() const {
    return seq_size_;
  }

  auto DiffBaseSeq() const {
    return diff_base_seq_;
  }

  auto DiffBaseSums() const {
    return diff_base_sums_;
  }

  auto MakeWrapper() const {
    return grammar::MakeDifferentialSLPWrapper(seq_size_,
                                               static_cast<const TSLP&>(*this),
                                               roots_,
                                               diff_base_seq_,
                                               span_sums_,
                                               diff_base_sums_,
                                               samples_,
                                               sample_roots_pos_,
                                               samples_pos_,
                                               samples_pos_rank_,
                                               samples_pos_select_);
  }

  void Compute(const sdsl::int_vector<>& da, uint32_t block_size);

  // Compute is split so the construct path can cache the bs-independent RePair
  // grammar and reuse it across all (bs,sf) cells:
  //  - BuildBaseGrammar: steps 1-2 (diff-encode + RePair into the TSLP base);
  //    sets seq_size_/diff_base_seq_; returns the top-level compact sequence.
  //  - RestoreBaseGrammar: install a cached base grammar + scalars (skips RePair).
  //  - FinishCompute: steps 3-6 (span sums + samples + field containers).
  std::vector<std::size_t> BuildBaseGrammar(const sdsl::int_vector<>& da);
  void RestoreBaseGrammar(const TSLP& base, std::size_t seq_size, uint64_t diff_base_seq);
  void FinishCompute(uint32_t block_size, const std::vector<std::size_t>& compact_seq);

  std::size_t serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, const std::string& name = "") const;

  void load(std::istream& in);

 private:
  TRoots roots_;
  TSpanSums span_sums_;
  TSamples samples_;
  TSampleRootsPos sample_roots_pos_;
  TBV samples_pos_;
  typename TBV::rank_1_type samples_pos_rank_;
  typename TBV::select_1_type samples_pos_select_;
  std::size_t seq_size_ = 0;
  uint64_t diff_base_seq_ = 0;
  uint64_t diff_base_sums_ = 0;
};

//~~~~~~~


template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
void DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::Compute(const sdsl::int_vector<>& da,
                                                                                       uint32_t block_size) {
  auto compact_seq = BuildBaseGrammar(da);
  FinishCompute(block_size, compact_seq);
}

// Steps 1-2: diff-encode the DA and RePair it into the TSLP base. bs-independent.
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
std::vector<std::size_t>
DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::BuildBaseGrammar(
    const sdsl::int_vector<>& da) {
  const auto n = da.size();
  seq_size_ = n;

  // 1. Compute diff-encoded DA: diff_da[0] = da[0] + diff_base_seq_; diff_da[i] = da[i]-da[i-1]+diff_base_seq_
  int min_diff = 0;
  for (std::size_t i = 1; i < n; ++i) {
    int diff = static_cast<int>(da[i]) - static_cast<int>(da[i - 1]);
    if (diff < min_diff)
      min_diff = diff;
  }
  diff_base_seq_ = (min_diff < 0) ? static_cast<uint64_t>(-min_diff) : 0;

  std::vector<int> diff_da(n);
  diff_da[0] = static_cast<int>(da[0]) + static_cast<int>(diff_base_seq_);
  for (std::size_t i = 1; i < n; ++i)
    diff_da[i] = static_cast<int>(da[i]) - static_cast<int>(da[i - 1]) + static_cast<int>(diff_base_seq_);

  // 2. Build non-CNF SLP from diff_da directly into TSLP base
  std::vector<std::size_t> compact_seq;
  {
    grammar::RePairEncoder<false> encoder;
    auto wrapper = grammar::BuildSLPWrapper(static_cast<TSLP&>(*this));
    auto report_c_seq = [&compact_seq](const auto& v) {
      compact_seq.emplace_back(v);
    };
    encoder.Encode(diff_da.begin(), diff_da.end(), wrapper, report_c_seq);
  }
  return compact_seq;
}

// Install a cached base grammar + scalars instead of re-running RePair.
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
void DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::RestoreBaseGrammar(
    const TSLP& base, std::size_t seq_size, uint64_t diff_base_seq) {
  static_cast<TSLP&>(*this) = base;
  seq_size_ = seq_size;
  diff_base_seq_ = diff_base_seq;
}

// Steps 3-6: span sums + samples + field containers, from the base + compact_seq.
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
void DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::FinishCompute(
    uint32_t block_size, const std::vector<std::size_t>& compact_seq) {
  const auto n = seq_size_;

  // Populate any cached SpanLength data that adapter TSLPs need (no-op by default).
  // Must run before ComputeSamplesOnCompactSequence, which reads SpanLength on roots.
  PopulateRootSpanLengths(static_cast<TSLP&>(*this), compact_seq);

  // 3. Compute span sums for non-terminals; ComputeSpanSums returns {min, max} of raw cover sums
  const auto& tslp = static_cast<const TSLP&>(*this);
  std::vector<uint64_t> raw_span_sums;
  auto report_span_sum = [&raw_span_sums](auto /*var*/, auto sum) {
    raw_span_sums.emplace_back(static_cast<uint64_t>(sum));
  };
  auto [min_ss, max_ss] = grammar::ComputeSpanSums(tslp, diff_base_seq_, report_span_sum);
  diff_base_sums_ = (min_ss < 0) ? static_cast<uint64_t>(-min_ss) : 0;

  // 4. Compute samples at block_size granularity
  const auto sigma = tslp.Sigma();
  const auto dbs_seq = diff_base_seq_;
  const auto dbs = diff_base_sums_;
  auto get_span_sum = [&tslp, &raw_span_sums, sigma, dbs_seq, dbs](auto var) -> int64_t {
    return tslp.IsTerminal(var) ? (static_cast<int64_t>(var) - static_cast<int64_t>(dbs_seq))
                                : (static_cast<int64_t>(raw_span_sums[var - sigma - 1]) - static_cast<int64_t>(dbs));
  };

  sdsl::bit_vector tmp_sample_pos(n, 0);
  std::vector<int64_t> tmp_samples;
  std::vector<std::size_t> tmp_roots_pos;
  auto report_sample = [&tmp_sample_pos, &tmp_samples, &tmp_roots_pos, n](auto pos, auto sum, auto root_idx) {
    if (pos < n)
      tmp_sample_pos[pos] = 1;
    tmp_samples.emplace_back(static_cast<int64_t>(sum));
    tmp_roots_pos.emplace_back(static_cast<std::size_t>(root_idx));
  };
  grammar::ComputeSamplesOnCompactSequence(compact_seq, tslp, get_span_sum, block_size, report_sample);

  // 5. Fill per-field int-vector containers. Built via a temporary sdsl::int_vector<>
  //    so compressed TSLP-side containers (enc_vector / dac_vector / vlc_vector) can be
  //    constructed from a range; for sdsl::int_vector<> the final ctor is an identity copy.
  auto fill_iv = [](auto& iv, const auto& vec) {
    using IV = std::decay_t<decltype(iv)>;
    sdsl::int_vector<> tmp(vec.size(), 0, 64);
    for (std::size_t i = 0; i < vec.size(); ++i)
      tmp[i] = static_cast<uint64_t>(vec[i]);
    sdsl::util::bit_compress(tmp);
    iv = IV(tmp);
  };

  fill_iv(roots_, compact_seq);
  fill_iv(span_sums_, raw_span_sums);
  fill_iv(samples_, tmp_samples);
  fill_iv(sample_roots_pos_, tmp_roots_pos);

  samples_pos_ = TBV(tmp_sample_pos);
  samples_pos_rank_ = typename TBV::rank_1_type(&samples_pos_);
  samples_pos_select_ = typename TBV::select_1_type(&samples_pos_);
}

//~~~~~~~


template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
std::size_t DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::serialize(
    std::ostream& out,
    sdsl::structure_tree_node* v,
    const std::string& name) const {
  std::size_t written = 0;
  written += TSLP::serialize(out);
  written += sdsl::serialize(roots_, out);
  written += sdsl::serialize(span_sums_, out);
  written += sdsl::serialize(samples_, out);
  written += sdsl::serialize(sample_roots_pos_, out);
  written += samples_pos_.serialize(out);
  // samples_pos_rank_ and samples_pos_select_ not serialized — rebuilt in load()
  written += sdsl::serialize(seq_size_, out);
  written += sdsl::serialize(diff_base_seq_, out);
  written += sdsl::serialize(diff_base_sums_, out);
  return written;
}

//~~~~~~~


template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
void DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::load(std::istream& in) {
  TSLP::load(in);
  sdsl::load(roots_, in);
  sdsl::load(span_sums_, in);
  sdsl::load(samples_, in);
  sdsl::load(sample_roots_pos_, in);
  samples_pos_.load(in);
  samples_pos_rank_ = typename TBV::rank_1_type(&samples_pos_);
  samples_pos_select_ = typename TBV::select_1_type(&samples_pos_);
  sdsl::load(seq_size_, in);
  sdsl::load(diff_base_seq_, in);
  sdsl::load(diff_base_sums_, in);
}

//~~~~~~~


// Per-field size breakdown for benchmark reporting.
template <typename TSLP,
          typename TRoots,
          typename TSpanSums,
          typename TSamples,
          typename TSampleRootsPos,
          typename TBV>
void collectSizes(SizeReport& out,
                  const DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& slp,
                  const std::string& prefix = "") {
  append(out, prefix + "base_slp", sdsl::size_in_bytes(static_cast<const TSLP&>(slp)));
  append(out, prefix + "roots", sdsl::size_in_bytes(slp.GetRoots()));
  append(out, prefix + "span_sums", sdsl::size_in_bytes(slp.GetSpanSums()));
  append(out, prefix + "samples", sdsl::size_in_bytes(slp.GetSamples()));
  append(out, prefix + "sample_roots_pos", sdsl::size_in_bytes(slp.GetSampleRootsPos()));
  append(out, prefix + "samples_pos", sdsl::size_in_bytes(slp.GetSamplesPos()));
}

//~~~~~~~


// ExpandSLP overload: more specific than the template in slp_tools.h, selected by overload resolution.
// Converts from exclusive-ep (ExpandSLP convention) to inclusive-ep (grammar::ExpandDifferentialSLP convention).
template <typename TSLP,
          typename TRoots,
          typename TSpanSums,
          typename TSamples,
          typename TSampleRootsPos,
          typename TBV,
          typename Report>
void ExpandSLP(const DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& slp,
               std::size_t bp,
               std::size_t ep,
               Report& report) {
  if (bp >= ep)
    return;
  auto wrapper = slp.MakeWrapper();
  grammar::ExpandDifferentialSLP(wrapper, bp, ep - 1, report);
}

//~~~~~~~


// On-disk cache of the bs-independent diff base grammar (TSLP base + the top-level
// compact sequence + the two diff scalars). Keyed by DiffGrammarCache<TSLP>'s type
// hash, so all (bs,sf) cells of a variant — and all container variants sharing the
// same TSLP — reuse one RePair pass.
template <typename TSLP>
struct DiffGrammarCache {
  using size_type = std::size_t;  // required by sdsl serialize/store_to_cache

  TSLP base;
  sdsl::int_vector<> compact_seq;
  uint64_t seq_size = 0;
  uint64_t diff_base_seq = 0;

  std::size_t serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const {
    auto* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    std::size_t b = 0;
    b += base.serialize(out, child, "base");
    b += compact_seq.serialize(out, child, "compact_seq");
    b += sdsl::write_member(seq_size, out, child, "seq_size");
    b += sdsl::write_member(diff_base_seq, out, child, "diff_base_seq");
    sdsl::structure_tree::add_size(child, b);
    return b;
  }
  void load(std::istream& in) {
    base.load(in);
    compact_seq.load(in);
    sdsl::read_member(seq_size, in);
    sdsl::read_member(diff_base_seq, in);
  }
};

// Install the diff base grammar into `work` (a DifferentialSLP/DifferentialLightSLP),
// loading the cached RePair grammar if present (bs-independent), else building +
// caching it. Returns the top-level compact sequence needed by FinishCompute.
template <typename TSLP, typename TDiff>
std::vector<std::size_t> LoadOrBuildDiffGrammar(TDiff& work, const sdsl::int_vector<>& da,
                                                Config& t_config, const std::string& grammar_key) {
  DiffGrammarCache<TSLP> gc;
  if (sdsl::cache_file_exists<DiffGrammarCache<TSLP>>(grammar_key, t_config)) {
    sdsl::load_from_cache(gc, grammar_key, t_config, true);
    work.RestoreBaseGrammar(gc.base, gc.seq_size, gc.diff_base_seq);
    return std::vector<std::size_t>(gc.compact_seq.begin(), gc.compact_seq.end());
  }
  auto compact_seq = work.BuildBaseGrammar(da);
  gc.base = static_cast<const TSLP&>(work);
  gc.compact_seq = sdsl::int_vector<>(compact_seq.size(), 0, 64);
  for (std::size_t i = 0; i < compact_seq.size(); ++i)
    gc.compact_seq[i] = compact_seq[i];
  sdsl::util::bit_compress(gc.compact_seq);
  gc.seq_size = work.SeqSize();
  gc.diff_base_seq = work.DiffBaseSeq();
  sdsl::store_to_cache(gc, grammar_key, t_config, true);
  return compact_seq;
}

template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples, typename TSampleRootsPos, typename TBV>
void construct(DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& t_dslp,
               Config& t_config,
               uint32_t block_size,
               const std::string& cache_key) {
  using namespace conf;
  using DSLP = DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>;

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  // Load-or-build the bs-independent RePair grammar (shared across container
  // variants), then finish the bs-dependent sampling. Convert into t_dslp to
  // bit-compress the base (per-field containers are already compressed).
  DSLP tmp;
  auto compact_seq = LoadOrBuildDiffGrammar<TSLP>(tmp, da, t_config, cache_key + "_grammar");
  tmp.FinishCompute(block_size, compact_seq);
  auto bit_compress = [](auto& v) {
    if constexpr (std::is_same_v<std::decay_t<decltype(v)>, sdsl::int_vector<>>)
      sdsl::util::bit_compress(v);
  };
  t_dslp = DSLP(tmp, bit_compress, bit_compress);

  sdsl::store_to_cache(t_dslp, cache_key, t_config, true);
}

}  // namespace dret
