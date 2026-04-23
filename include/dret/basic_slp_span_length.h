//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/23/26.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <set>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <grammar/slp.h>
#include <grammar/utility.h>

namespace dret {

// Default customization point: for any TSLP that already provides SpanLength
// natively (e.g. grammar::SLP<>), do nothing. Adapters that need to populate a
// length cache from compact_seq (DifferentialSLP's roots_) overload this.
template <typename TSLP, typename RootsSeq>
void PopulateRootSpanLengths(TSLP& /*slp*/, const RootsSeq& /*compact_seq*/) {}

//~~~~~~~


// On-the-fly SpanLength: inherits from a BasicSLP-like type and computes
// SpanLength by recursive descent to terminals. No stored lengths_; cost per
// call is O(span).
template <typename TBasicSLP = grammar::BasicSLP<>>
class BasicSLPOnTheFlySpanLength : public TBasicSLP {
 public:
  using VariableType = typename TBasicSLP::VariableType;

  using TBasicSLP::TBasicSLP;

  BasicSLPOnTheFlySpanLength() = default;

  // Mirrors grammar::SLP's 2-action copy-ctor signature so DifferentialSLP's
  // templated copy-ctor (rules + lengths + int-containers actions) keeps
  // compiling. The lengths action is a no-op here.
  template <typename TOther,
            typename ActionVars = grammar::NoAction,
            typename ActionLengths = grammar::NoAction>
  BasicSLPOnTheFlySpanLength(const BasicSLPOnTheFlySpanLength<TOther>& other,
                             ActionVars&& action_vars = grammar::NoAction(),
                             ActionLengths&& /*action_lengths*/ = grammar::NoAction())
      : TBasicSLP(static_cast<const TOther&>(other), std::forward<ActionVars>(action_vars)) {}

  std::size_t SpanLength(VariableType var) const {
    if (this->IsTerminal(var))
      return 1;
    auto children = (*this)[var];
    return SpanLength(children.first) + SpanLength(children.second);
  }
};

//~~~~~~~


// Subset-lengths SpanLength: stores lengths only for non-terminals that appear
// in compact_seq (DifferentialSLP's roots). Hit: O(1) via bit_vector + rank.
// Miss: O(span) fallback descent — used for interior left children during
// ExpandDifferentialSLP's recursive descent.
template <typename TBasicSLP = grammar::BasicSLP<>,
          typename TIntContainer = sdsl::int_vector<>,
          typename TBV = sdsl::sd_vector<>>
class BasicSLPCachedRootSpanLengths : public TBasicSLP {
 public:
  using VariableType = typename TBasicSLP::VariableType;

  using TBasicSLP::TBasicSLP;

  BasicSLPCachedRootSpanLengths() = default;

  template <typename TOtherBasicSLP,
            typename TOtherIntContainer,
            typename TOtherBV,
            typename ActionVars = grammar::NoAction,
            typename ActionLengths = grammar::NoAction,
            typename ActionInt = grammar::NoAction>
  BasicSLPCachedRootSpanLengths(
      const BasicSLPCachedRootSpanLengths<TOtherBasicSLP, TOtherIntContainer, TOtherBV>& other,
      ActionVars&& action_vars = grammar::NoAction(),
      ActionLengths&& /*action_lengths*/ = grammar::NoAction(),
      ActionInt&& action_int = grammar::NoAction())
      : TBasicSLP(static_cast<const TOtherBasicSLP&>(other), std::forward<ActionVars>(action_vars)) {
    grammar::Construct(cached_lengths_, other.GetCachedLengths());
    action_int(cached_lengths_);
    grammar::Construct(cached_mask_, other.GetCachedMask());
    cached_rank_ = typename TBV::rank_1_type(&cached_mask_);
  }

  template <typename RootsSeq>
  void PopulateCache(const RootsSeq& compact_seq) {
    const auto sigma = this->Sigma();
    const auto n_nt = this->Variables() - sigma;

    std::set<std::uint64_t> uniq;
    for (const auto& v : compact_seq) {
      const auto var = static_cast<std::uint64_t>(v);
      if (!this->IsTerminal(static_cast<VariableType>(var)))
        uniq.insert(var);
    }

    sdsl::bit_vector mask(n_nt, 0);
    for (auto var : uniq) {
      const auto idx = var - static_cast<std::uint64_t>(sigma) - 1;
      mask[idx] = 1;
    }

    std::vector<std::uint64_t> lengths;
    lengths.reserve(uniq.size());
    for (auto var : uniq) {
      lengths.emplace_back(computeSpanLength(static_cast<VariableType>(var)));
    }

    cached_mask_ = TBV(mask);
    cached_rank_ = typename TBV::rank_1_type(&cached_mask_);

    cached_lengths_ = TIntContainer(lengths.size(), 0, 64);
    for (std::size_t i = 0; i < lengths.size(); ++i)
      cached_lengths_[i] = lengths[i];
    sdsl::util::bit_compress(cached_lengths_);
  }

  std::size_t SpanLength(VariableType var) const {
    if (this->IsTerminal(var))
      return 1;
    const auto idx = static_cast<std::uint64_t>(var) - static_cast<std::uint64_t>(this->Sigma()) - 1;
    if (idx < cached_mask_.size() && cached_mask_[idx])
      return static_cast<std::size_t>(cached_lengths_[cached_rank_(idx)]);
    return computeSpanLength(var);
  }

  const auto& GetCachedLengths() const {
    return cached_lengths_;
  }

  const auto& GetCachedMask() const {
    return cached_mask_;
  }

  std::size_t serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, const std::string& name = "") const {
    std::size_t written = 0;
    written += TBasicSLP::serialize(out);
    written += sdsl::serialize(cached_lengths_, out);
    written += cached_mask_.serialize(out);
    return written;
  }

  void load(std::istream& in) {
    TBasicSLP::load(in);
    sdsl::load(cached_lengths_, in);
    cached_mask_.load(in);
    cached_rank_ = typename TBV::rank_1_type(&cached_mask_);
  }

 private:
  std::size_t computeSpanLength(VariableType var) const {
    if (this->IsTerminal(var))
      return 1;
    auto children = (*this)[var];
    return computeSpanLength(children.first) + computeSpanLength(children.second);
  }

  TIntContainer cached_lengths_;
  TBV cached_mask_;
  typename TBV::rank_1_type cached_rank_;
};

//~~~~~~~


template <typename TBasicSLP, typename TIntContainer, typename TBV, typename RootsSeq>
void PopulateRootSpanLengths(BasicSLPCachedRootSpanLengths<TBasicSLP, TIntContainer, TBV>& slp,
                             const RootsSeq& compact_seq) {
  slp.PopulateCache(compact_seq);
}

}  // namespace dret
