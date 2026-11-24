//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/9/20.
//

#ifndef DRET_INCLUDE_DRET_DIFF_SLP_H_
#define DRET_INCLUDE_DRET_DIFF_SLP_H_

#include <sdsl/bit_vectors.hpp>

#include <grammar/slp.h>

namespace grammar {

//template<typename SLP, typename CompactSeq, typename GetCoverSum, typename ReportAccumulatedSum>
//void ComputeSampleSumsOnCompactSequence(const CompactSeq &_compact_seq,
//                                        const SLP &_slp,
//                                        const GetCoverSum &_get_cover_sum,
//                                        const ReportAccumulatedSum &_report) {
//  int64_t sum = 0;
//  std::size_t pos = 0;
//
//  for (std::size_t i = 0; i < _compact_seq.size(); ++i) {
//    pos += _slp.SpanLength(_compact_seq[i]);
//    sum += _get_cover_sum(_compact_seq[i]);
//
//    _report(pos, sum);
//  }
//}

template<typename SLP, typename CompactSeq, typename DifferentialBase, typename ReportCoverSum, typename ReportAccumulatedSum>
void ComputeCoverSumsAndSampleSum(const SLP &_slp,
                                  const CompactSeq &_compact_seq,
                                  const DifferentialBase &_diff_base,
                                  ReportCoverSum &_report_cover_sum,
                                  ReportAccumulatedSum &_report_acc_sum) {
  auto rules = _slp.GetRules();
  std::vector<int64_t> cover_sums(rules.size());

  auto sigma = _slp.Sigma();
  uint32_t diff_base = 0; // It is not considered for partial sums in each variable.
  auto get_cover_sum = [&_slp, &diff_base, &cover_sums, &sigma](auto _var) {
    return _slp.IsTerminal(_var) ? (_var - diff_base) : (cover_sums[_var - sigma - 1] - _slp.SpanLength(_var) * diff_base);
  };

  // Compute cover sum for each non-terminal
  for (std::size_t i = 0, j = sigma + 1; i < rules.size() / 2; ++i, ++j) {
    auto children = _slp[j];
    cover_sums[i] = get_cover_sum(children.first) + get_cover_sum(children.second);
    _report_cover_sum(i, cover_sums[i]);
  }

  // Compute accumulated sums on compact sequence
  uint32_t sum = 0;
  std::size_t pos = 0;
  diff_base = _diff_base; // Now, it is considered for accumulated sample sums.
  for (std::size_t i = 0; i < _compact_seq.size(); ++i) {
    _report_acc_sum(pos, sum);

    pos += _slp.SpanLength(_compact_seq[i]);
    sum += get_cover_sum(_compact_seq[i]);
  }
}

template<typename SLP = grammar::SLP<>,
    typename Roots = std::vector<std::size_t>,
    typename SpanSums = std::vector<uint32_t>,
    typename SampleSums = std::vector<uint32_t>,
    typename BitVector = sdsl::sd_vector<>,
    typename BitVectorRank = typename BitVector::rank_1_type,
    typename BitVectorSelect = typename BitVector::select_1_type>
class DifferentialSLP : public SLP {
 public:
  using size_type = std::size_t;

  DifferentialSLP() = default;

  template<typename OriginalSLP, typename CompactSeq, typename SpanSumsAction = NoAction, typename SampleSumsAction = NoAction, typename RootsAction = NoAction>
  void Compute(std::size_t _seq_size,
               const OriginalSLP &_slp,
               const CompactSeq &_compact_seq,
               uint64_t _differential_base,
               SpanSumsAction &&_span_sums_action = NoAction(),
               SampleSumsAction &&_sample_sums_action = NoAction(),
               RootsAction &&_roots_action = NoAction()) {

    differential_base_ = _differential_base;

    std::vector<typename SpanSums::value_type> span_sums;
    span_sums.reserve(_slp.GetRules().size());
    auto report_cover_sum = [&span_sums](auto _rule, auto _sum) {
      span_sums.emplace_back(_sum);
    };

    sdsl::bit_vector tmp_sample_pos(_seq_size + 1, 0);
    std::vector<typename SampleSums::value_type> sample_sums;
    auto report_acc_sum = [&tmp_sample_pos, &sample_sums](auto _pos, auto _sum) {
      tmp_sample_pos[_pos] = 1;
      sample_sums.emplace_back(_sum);
    };

    ComputeCoverSumsAndSampleSum(_slp, _compact_seq, _differential_base, report_cover_sum, report_acc_sum);
    tmp_sample_pos[tmp_sample_pos.size() - 1] = 1;

    OriginalSLP::operator=(_slp);

//    span_sums_ = SpanSums(span_sums);
    Construct(span_sums_, span_sums);
    _span_sums_action(span_sums_);
    span_sums.clear();

//    sample_sums_ = SampleSums(sample_sums);
    Construct(sample_sums_, sample_sums);
    _sample_sums_action(sample_sums_);
    sample_sums.clear();

    sample_pos_ = BitVector(tmp_sample_pos);
    sample_pos_rank_ = BitVectorRank(&sample_pos_);
    sample_pos_select_ = BitVectorSelect(&sample_pos_);

//    roots_ = Roots(_compact_seq);
    Construct(roots_, _compact_seq);
    _roots_action(roots_);

  }

  auto DifferentialBase() const {
    return differential_base_;
  }

  auto SpanSum(std::size_t _var) const {
    return SLP::IsTerminal(_var) ? _var : span_sums_[_var - SLP::sigma_ - 1];
  }

  auto Root(std::size_t _i) const {
    assert(_i < roots_.size());
    return roots_[_i];
  }

  std::size_t Sample(std::size_t _pos) const {
    return sample_pos_rank_(_pos + 1);
  }

  auto SamplePosition(std::size_t _sample) const {
    return sample_pos_select_(_sample);
  }

  auto SampleFirstRoot(std::size_t _sample) const {
    return _sample - 1;
  }

  auto SampleSum(std::size_t _sample) const {
    assert(0 < _sample);
    return sample_sums_[_sample - 1];
  }

  std::size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, const std::string &name = "") const {
    std::size_t written_bytes = 0;
    written_bytes += SLP::serialize(out);

    written_bytes += sdsl::serialize(roots_, out);

    written_bytes += sdsl::serialize(differential_base_, out);
    written_bytes += sdsl::serialize(span_sums_, out);

    written_bytes += sample_pos_.serialize(out);
    written_bytes += sample_pos_rank_.serialize(out);
    written_bytes += sample_pos_select_.serialize(out);
    written_bytes += sdsl::serialize(sample_sums_, out);

    return written_bytes;
  }

  void load(std::istream &in) {
    SLP::load(in);

    sdsl::load(roots_, in);

    sdsl::load(differential_base_, in);
    sdsl::load(span_sums_, in);

    sample_pos_.load(in);
//    sample_pos_rank_.load(in);
    sample_pos_rank_ = BitVectorRank(&sample_pos_);
//    sample_pos_select_.load(in);
    sample_pos_select_ = BitVectorSelect(&sample_pos_);
    sdsl::load(sample_sums_, in);
  }

 private:
  Roots roots_; // Compact sequence. It is a sequence of variables (terminals/non-terminals), i.e. a forest.

  uint64_t differential_base_;
  SpanSums span_sums_; // For each non-terminal.

  BitVector sample_pos_; // Marks sampled positions in the original sequence.
  BitVectorRank sample_pos_rank_;
  BitVectorSelect sample_pos_select_;
  SampleSums sample_sums_;
};

template<typename DiffSLP, typename Report>
void ExpandSLPFromLeft(const DiffSLP &_slp, std::size_t _var, std::size_t &_length, const Report &_report) {
  assert(0 < _length);

  if (_slp.IsTerminal(_var)) {
    _report(_var);
    --_length;
    return;
  }

  const auto &children = _slp[_var];
  ExpandSLPFromLeft(_slp, children.first, _length, _report);

  if (0 < _length) {
    ExpandSLPFromLeft(_slp, children.second, _length, _report);
  }
}

template<typename DiffSLP, typename Report, typename Skip>
void ExpandSLPFromLeft(const DiffSLP &_slp,
                       std::size_t _var,
                       std::size_t &_sp,
                       std::size_t &_length,
                       const Report &_report,
                       const Skip &_skip) {
  assert(0 < _sp);
  assert(0 < _length);

//  if (_slp.IsTerminal(_var)) {
//    _skip(_var);
//    --_sp;
//    return;
//  }
//
//  const auto &children = _slp[_var];
//  ExpandSLPFromLeft(_slp, children.first, _sp, _length, _report, _skip);

  const auto &children = _slp[_var];
  auto left_child_len = _slp.SpanLength(children.first);
  if (_sp < left_child_len) {
    ExpandSLPFromLeft(_slp, children.first, _sp, _length, _report, _skip);
  } else {
    _skip(children.first);
    _sp -= left_child_len;
  }

  if (0 < _sp) {
    ExpandSLPFromLeft(_slp, children.second, _sp, _length, _report, _skip);
  } else if (0 < _length) {
    ExpandSLPFromLeft(_slp, children.second, _length, _report);
  }
}

template<typename DiffSLP, typename Report, typename Skip>
void ExpandSLPFromFront(const DiffSLP &_slp,
                        std::size_t _idx_root,
                        std::size_t _sp,
                        std::size_t _length,
                        const Report &_report,
                        const Skip &_skip) {
  assert(0 < _length);

  std::size_t root;
  std::size_t span_length;
  while (root = _slp.Root(_idx_root), (span_length = _slp.SpanLength(root)) <= _sp) {
    _skip(root);

    _sp -= span_length;
    ++_idx_root;
  }

  if (0 < _sp) {
    ExpandSLPFromLeft(_slp, root, _sp, _length, _report, _skip);
    ++_idx_root;
  }

  while (0 < _length) {
    ExpandSLPFromLeft(_slp, _slp.Root(_idx_root), _length, _report);
    ++_idx_root;
  }
}

template<typename DiffSLP, typename Report>
void ExpandDifferentialSLP(const DiffSLP &_slp, std::size_t _sp, std::size_t _ep, const Report &_report) {
  assert(_sp <= _ep);

  auto sample = _slp.Sample(_sp);
  auto idx_root = _slp.SampleFirstRoot(sample);
  auto pos = _slp.SamplePosition(sample);
  auto sp = _sp - pos;
  auto len = _ep - _sp + 1;

  auto sum = _slp.SampleSum(sample);
  auto diff_base = _slp.DifferentialBase();
  auto report = [&_report, &sum, &diff_base](const auto &_terminal) {
    sum += _terminal - diff_base;
    _report(sum);
  };

  auto skip = [&sum, &_slp, &diff_base](const auto _var) {
    sum += _slp.SpanSum(_var) - _slp.SpanLength(_var) * diff_base;
  };

  ExpandSLPFromFront(_slp, idx_root, sp, len, report, skip);
}

} // namespace grammar
#endif //DRET_INCLUDE_DRET_DIFF_SLP_H_
