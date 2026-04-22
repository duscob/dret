//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/22/26.
//

#pragma once

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <grammar/differential_slp.h>
#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "config.h"

namespace dret {

template <typename TSLP          = grammar::SLP<>,
          typename TSampledSLP   = grammar::SampledSLP<>,
          typename TIntContainer = sdsl::int_vector<>,
          typename TBV           = sdsl::sd_vector<>>
class DifferentialLightSLP : public TSLP, public TSampledSLP {
 public:
  using size_type = std::size_t;

  auto MakeWrapper() const {
    return grammar::MakeDifferentialSLPWrapper(
        seq_size_,
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

  void Compute(const sdsl::int_vector<>& da,
               uint32_t block_size,
               float storing_factor,
               grammar::Chunks<>& cslp_docs_out);

  std::size_t serialize(std::ostream& out,
                        sdsl::structure_tree_node* v = nullptr,
                        const std::string& name = "") const;

  void load(std::istream& in);

 private:
  TIntContainer               roots_;
  TIntContainer               span_sums_;
  TIntContainer               samples_;
  TIntContainer               sample_roots_pos_;
  TBV                         samples_pos_;
  typename TBV::rank_1_type   samples_pos_rank_;
  typename TBV::select_1_type samples_pos_select_;
  std::size_t                 seq_size_       = 0;
  uint64_t                    diff_base_seq_  = 0;
  uint64_t                    diff_base_sums_ = 0;
};

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV>
void DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>::Compute(
    const sdsl::int_vector<>& da,
    uint32_t block_size,
    float storing_factor,
    grammar::Chunks<>& cslp_docs_out) {
  const auto n = da.size();
  seq_size_ = n;

  // 1. Compute diff-encoded DA: diff_da[0] = da[0] + diff_base_seq_; diff_da[i] = da[i]-da[i-1]+diff_base_seq_
  int min_diff = 0;
  for (std::size_t i = 1; i < n; ++i) {
    int diff = static_cast<int>(da[i]) - static_cast<int>(da[i - 1]);
    if (diff < min_diff) min_diff = diff;
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
    auto report_c_seq = [&compact_seq](const auto& v) { compact_seq.emplace_back(v); };
    encoder.Encode(diff_da.begin(), diff_da.end(), wrapper, report_c_seq);
  }

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
    return tslp.IsTerminal(var)
        ? (static_cast<int64_t>(var) - static_cast<int64_t>(dbs_seq))
        : (static_cast<int64_t>(raw_span_sums[var - sigma - 1]) - static_cast<int64_t>(dbs));
  };

  sdsl::bit_vector tmp_sample_pos(n, 0);
  std::vector<int64_t> tmp_samples;
  std::vector<std::size_t> tmp_roots_pos;
  auto report_sample = [&tmp_sample_pos, &tmp_samples, &tmp_roots_pos, n](auto pos, auto sum, auto root_idx) {
    if (pos < n) tmp_sample_pos[pos] = 1;
    tmp_samples.emplace_back(static_cast<int64_t>(sum));
    tmp_roots_pos.emplace_back(static_cast<std::size_t>(root_idx));
  };
  grammar::ComputeSamplesOnCompactSequence(compact_seq, tslp, get_span_sum, block_size, report_sample);

  // 5. Fill bit-compressed TIntContainer fields
  auto fill_iv = [](TIntContainer& iv, const auto& vec) {
    iv = TIntContainer(vec.size(), 0, 64);
    for (std::size_t i = 0; i < vec.size(); ++i)
      iv[i] = static_cast<uint64_t>(vec[i]);
    sdsl::util::bit_compress(iv);
  };

  fill_iv(roots_, compact_seq);
  fill_iv(span_sums_, raw_span_sums);
  fill_iv(samples_, tmp_samples);  // accumulated sums >= 0 at sample points (DA values)
  fill_iv(sample_roots_pos_, tmp_roots_pos);

  samples_pos_ = TBV(tmp_sample_pos);
  samples_pos_rank_ = typename TBV::rank_1_type(&samples_pos_);
  samples_pos_select_ = typename TBV::select_1_type(&samples_pos_);

  // 6. Build SampledSLP base from original DA using CombinedSLP (for computeCover)
  {
    std::vector<int> da_vec(n);
    for (std::size_t i = 0; i < n; ++i) da_vec[i] = static_cast<int>(da[i]);

    grammar::SLP<> slp_cnf;
    {
      grammar::RePairEncoder<true> encoder;
      auto wrapper = grammar::BuildSLPWrapper(slp_cnf);
      encoder.Encode(da_vec.begin(), da_vec.end(), wrapper);
    }

    grammar::CombinedSLP<> cslp(slp_cnf);
    grammar::AddSet add_set(cslp_docs_out);
    cslp.Compute(block_size, add_set, add_set,
                 grammar::MustBeSampled<grammar::Chunks<>>(
                     grammar::AreChildrenTooBig(cslp_docs_out, storing_factor)));

    // Assign TSampledSLP base; rank/select pointers temporarily dangling — fixed on load from cache
    static_cast<TSampledSLP&>(*this) = static_cast<const TSampledSLP&>(cslp);
  }
}

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV>
std::size_t DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>::serialize(
    std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const {
  std::size_t written = 0;
  written += TSLP::serialize(out);
  written += TSampledSLP::serialize(out);
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


template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV>
void DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>::load(std::istream& in) {
  TSLP::load(in);
  TSampledSLP::load(in);
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


// ExpandSLP overload: more specific than the template in slp_tools.h, selected by overload resolution.
// Converts from exclusive-ep (ExpandSLP convention) to inclusive-ep (grammar::ExpandDifferentialSLP convention).
template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV, typename Report>
void ExpandSLP(const DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>& slp,
               std::size_t bp, std::size_t ep, Report& report) {
  if (bp >= ep) return;
  auto wrapper = slp.MakeWrapper();
  grammar::ExpandDifferentialSLP(wrapper, bp, ep - 1, report);
}

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV>
void construct(DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>& t_dslp,
               Config& t_config, uint32_t block_size, float storing_factor) {
  using namespace conf;

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  grammar::Chunks<> cslp_docs;
  t_dslp.Compute(da, block_size, storing_factor, cslp_docs);

  const std::string key_prefix = std::format("{}-{}_", block_size, storing_factor);
  const std::string key_docs = key_prefix + t_config.keys[kDGCDA][kDocs].get<std::string>();

  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);
  auto bit_compress = [](sdsl::int_vector<>& v) { sdsl::util::bit_compress(v); };
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);

  const std::string key_slp = key_prefix + t_config.keys[kDGCDA][kSLP].get<std::string>();
  sdsl::store_to_cache(t_dslp, key_slp, t_config, true);
}

}  // namespace dret
