//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// Shared RMQ-based machinery used by both document-frequency and document-listing
// indexes (doc_freq_index_rmq.h and doc_list_index_rmq.h).
//

#pragma once

#include <cstddef>
#include <stack>
#include <utility>

namespace dret {

enum class OccurrenceSide { LEFTMOST, RIGHTMOST };

// TODO Specialize to avoid conditionals in LEFTMOST or RIGHTMOST cases.
template <OccurrenceSide _side,
          typename SetFirstRanges,
          typename RMQ,
          typename GetValue,
          typename IsReported,
          typename Report>
void GetExtremeOccurrencesRMQ(std::size_t _sp,
                              std::size_t _ep,
                              const SetFirstRanges& _set_first_ranges,
                              const RMQ& _rmq,
                              const GetValue& _get_value,
                              const IsReported& _is_reported,
                              Report& _report) {
  using Range = std::pair<std::size_t, std::size_t>;
  std::stack<Range> stack;

  _set_first_ranges(_sp, _ep, stack, _side);

  while (!stack.empty()) {
    auto values = stack.top();
    auto& rsp = values.first;
    auto& rep = values.second;

    stack.pop();

    if (rsp <= rep) {
      std::size_t idx = _rmq(rsp, rep);

      auto value = _get_value(idx, _side, _sp, _ep);

      if (!_is_reported(idx, value, _side)) {
        _report(idx, value, _side, _sp, _ep);

        // For leftmost occurrences, we must search in the right range after the left range.
        if (_side == OccurrenceSide::LEFTMOST) stack.emplace(idx + 1, rep);

        // Guard against unsigned underflow when idx == 0 (the original SADA-on-SA
        // setup avoided it via the `\0` sentinel at SA[0], but the ILCP run-coord
        // recursion can land at idx == 0 legitimately).
        if (idx > 0) stack.emplace(rsp, idx - 1);

        // For rightmost occurrences, we must search in the right range before in the left range.
        if (_side == OccurrenceSide::RIGHTMOST) stack.emplace(idx + 1, rep);
      }
    }
  }
}

template <typename LeftRMQ, typename RightRMQ>
struct RMQAlgoCore {
 public:
  LeftRMQ left_rmq;
  RightRMQ right_rmq;

  // Set initial ranges in stack
  using Stack = std::stack<std::pair<std::size_t, std::size_t>>;
  virtual void operator()(std::size_t _sp, std::size_t _ep, Stack& _stack, OccurrenceSide _side) const = 0;
};

template <typename LeftRMQ, typename RightRMQ>
struct RMQAlgoCoreSada : public RMQAlgoCore<LeftRMQ, RightRMQ> {
 public:
  using typename RMQAlgoCore<LeftRMQ, RightRMQ>::Stack;
  void operator()(std::size_t _sp, std::size_t _ep, Stack& _stack, OccurrenceSide) const override {
    _stack.emplace(_sp, _ep);
  }
};

template <typename LeftRMQ,
          typename RightRMQ,
          typename BitVector,
          typename BitVectorRank = typename BitVector::rank_1_type,
          typename BitVectorSelect = typename BitVector::select_1_type>
struct RMQAlgoCoreILCP : public RMQAlgoCore<LeftRMQ, RightRMQ> {
 public:
  BitVector run_heads[2];
  BitVectorRank run_heads_rank[2];
  BitVectorSelect run_heads_select[2];

  using typename RMQAlgoCore<LeftRMQ, RightRMQ>::Stack;
  void operator()(std::size_t _sp, std::size_t _ep, Stack& _stack, OccurrenceSide _side) const override {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

    _stack.emplace(run_heads_rank[direction](_sp + 1) - 1, run_heads_rank[direction](_ep + 1) - 1);
  }
};

template <typename LeftRMQ,
          typename RightRMQ,
          typename BitVector,
          typename BitVectorRank = typename BitVector::rank_1_type,
          typename BitVectorSelect = typename BitVector::select_1_type>
struct RMQAlgoCoreCILCP : public RMQAlgoCoreILCP<LeftRMQ, RightRMQ, BitVector, BitVectorRank, BitVectorSelect> {
 public:
  using typename RMQAlgoCore<LeftRMQ, RightRMQ>::Stack;
  void operator()(std::size_t _sp, std::size_t _ep, Stack& _stack, OccurrenceSide _side) const override {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

    auto p = std::make_pair(this->run_heads_rank[direction](_sp + 1) - 1,
                            this->run_heads_rank[direction](_ep + 1) - 1);

    if (p.first < p.second) {
      if (_side == dret::OccurrenceSide::LEFTMOST) {
        _stack.emplace(p.second, p.second);
        _stack.emplace(p.first, p.second - 1);
      } else if (_side == dret::OccurrenceSide::RIGHTMOST) {
        _stack.emplace(p.first, p.first);
        _stack.emplace(p.first + 1, p.second);
      }
    } else {
      _stack.emplace(p);
    }
  }
};

}  // namespace dret
