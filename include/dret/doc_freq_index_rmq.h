//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/6/20.
//

#ifndef DRET_DOC_FREQ_INDEX_RMQ_H_
#define DRET_DOC_FREQ_INDEX_RMQ_H_

#include <cstddef>
#include <utility>
#include <stack>
#include <vector>
#include <algorithm>

#include "doc_freq_index.h"

namespace dret {

enum class OccurrenceSide { LEFTMOST, RIGHTMOST };

// TODO Specialize to avoid conditionals in LEFTMOST or RIGHTMOST cases.
template<OccurrenceSide _side, typename SetFirstRanges, typename RMQ, typename GetValue, typename IsReported, typename Report>
void GetExtremeOccurrencesRMQ(std::size_t _sp,
                              std::size_t _ep,
                              const SetFirstRanges &_set_first_ranges,
                              const RMQ &_rmq,
                              const GetValue &_get_value,
                              const IsReported &_is_reported,
                              Report &_report) {
  using Range = std::pair<std::size_t, std::size_t>;
  std::stack<Range> stack;

  _set_first_ranges(_sp, _ep, stack, _side);

  while (!stack.empty()) {
    auto values = stack.top();
    auto &rsp = values.first;
    auto &rep = values.second;

    stack.pop();

    if (rsp <= rep) {
      std::size_t idx = _rmq(rsp, rep);

      auto value = _get_value(idx, _side, _sp, _ep);

      if (!_is_reported(idx, value, _side)) {
        _report(idx, value, _side, _sp, _ep);

        // For leftmost occurrences, we must search in the right range after the left range.
        if (_side == OccurrenceSide::LEFTMOST)
          stack.emplace(idx + 1, rep);

        stack.emplace(rsp, idx - 1); // idx != 0, since `\0` is appended to string

        // For rightmost occurrences, we must search in the right range before in the left range.
        if (_side == OccurrenceSide::RIGHTMOST)
          stack.emplace(idx + 1, rep);
      }
    }
  }
}

template<typename LeftRMQ, typename RightRMQ, typename SetFirstRange, typename GetSuffixDocValue, typename ReportSuffixDocValue, typename IsReported, typename MarkDoc>
auto ComputeSuffixesRMQ(std::size_t _sp,
                        std::size_t _ep,
                        const LeftRMQ &_left_rmq,
                        const RightRMQ &_right_rmq,
                        const SetFirstRange &_set_first_range,
                        const GetSuffixDocValue &_get_suffix_doc_values,
                        const ReportSuffixDocValue &_report_suffix_doc_value,
                        const IsReported &_is_reported,
                        MarkDoc _mark_doc) {
  std::vector<std::size_t> suffixes;
  auto add_suffix = [&suffixes](auto _suffix) {
    suffixes.emplace_back(_suffix);
  };

  auto report_suffix =
      [&_report_suffix_doc_value, &add_suffix, &_mark_doc](
          auto _i, const auto &_suffix_doc, OccurrenceSide _side, auto _sp, auto _ep) {
        _report_suffix_doc_value(_i, _suffix_doc, _side, add_suffix, _mark_doc, _sp, _ep);
      };

  GetExtremeOccurrencesRMQ<OccurrenceSide::LEFTMOST>(
      _sp, _ep, _set_first_range, _left_rmq, _get_suffix_doc_values, _is_reported, report_suffix);

  GetExtremeOccurrencesRMQ<OccurrenceSide::RIGHTMOST>(
      _sp, _ep, _set_first_range, _right_rmq, _get_suffix_doc_values, _is_reported, report_suffix);

  std::sort(suffixes.begin(), suffixes.end());

  return suffixes;
}

template<typename LeftRMQ, typename RightRMQ>
struct RMQAlgoCore {
 public:
//  std::size_t doc_cnt;

  LeftRMQ left_rmq;
  RightRMQ right_rmq;

  // Set initial ranges in stack
  using Stack = std::stack<std::pair<std::size_t, std::size_t>>;
  virtual void operator()(std::size_t _sp, std::size_t _ep, Stack &_stack, OccurrenceSide _side) const = 0;
};

template<typename LeftRMQ, typename RightRMQ>
struct RMQAlgoCoreSada : public RMQAlgoCore<LeftRMQ, RightRMQ> {
 public:
  // Set initial ranges
  using typename RMQAlgoCore<LeftRMQ, RightRMQ>::Stack;
  void operator()(std::size_t _sp, std::size_t _ep, Stack &_stack, OccurrenceSide) const override {
    _stack.emplace(_sp, _ep);
  }
};

template<typename LeftRMQ, typename RightRMQ, typename BitVector, typename BitVectorRank = typename BitVector::rank_1_type, typename BitVectorSelect = typename BitVector::select_1_type>
struct RMQAlgoCoreILCP : public RMQAlgoCore<LeftRMQ, RightRMQ> {
 public:
  BitVector run_heads[2];
  BitVectorRank run_heads_rank[2];
  BitVectorSelect run_heads_select[2];

  // Set initial ranges
  using typename RMQAlgoCore<LeftRMQ, RightRMQ>::Stack;
  void operator()(std::size_t _sp, std::size_t _ep, Stack &_stack, OccurrenceSide _side) const override {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

    _stack.emplace(run_heads_rank[direction](_sp + 1) - 1, run_heads_rank[direction](_ep + 1) - 1);
  }
};

template<typename LeftRMQ, typename RightRMQ, typename BitVector, typename BitVectorRank = typename BitVector::rank_1_type, typename BitVectorSelect = typename BitVector::select_1_type>
struct RMQAlgoCoreCILCP : public RMQAlgoCoreILCP<LeftRMQ, RightRMQ, BitVector, BitVectorRank, BitVectorSelect> {
 public:
  // Set initial ranges
  using typename RMQAlgoCore<LeftRMQ, RightRMQ>::Stack;
  void operator()(std::size_t _sp, std::size_t _ep, Stack &_stack, OccurrenceSide _side) const override {
    auto direction = _side == dret::OccurrenceSide::LEFTMOST ? 0 : 1;

    auto p = std::make_pair(this->run_heads_rank[direction](_sp + 1) - 1, this->run_heads_rank[direction](_ep + 1) - 1);

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

class GetValuesRMQFunctor {
 public:
  virtual std::pair<std::size_t, std::size_t> operator()(std::size_t _i,
                                                         const OccurrenceSide &_side,
                                                         std::size_t _sp,
                                                         std::size_t _ep) const = 0;
};

template<typename GetValue>
class GetValuesRMQSadaFunctor : public GetValuesRMQFunctor {
 public:
  GetValuesRMQSadaFunctor(const GetValue *_get_value) : get_value_{_get_value} {}

  std::pair<std::size_t, std::size_t> operator()(std::size_t _i,
                                                 const OccurrenceSide &_side,
                                                 std::size_t _sp,
                                                 std::size_t _ep) const override {
    return (*get_value_)(_i);
  }

 private:
  const GetValue *get_value_;
};

class RMQReporter {
 public:
  virtual void operator()(std::size_t _i,
                          const std::pair<std::size_t, std::size_t> &_suffix_doc,
                          const OccurrenceSide &_side,
                          const std::function<void(std::size_t)> &_report,
                          std::size_t _sp,
                          std::size_t _ep) const = 0;
};

template<typename Mark>
class RMQSadaReporter : public RMQReporter {
 public:
  RMQSadaReporter(Mark *_mark) : mark_{_mark} {}

  void operator()(std::size_t _i,
                  const std::pair<std::size_t, std::size_t> &_suffix_doc,
                  const OccurrenceSide &_side,
                  const std::function<void(std::size_t)> &_report,
                  std::size_t _sp,
                  std::size_t _ep) const override {
    _report(_suffix_doc.first);
    (*mark_)(_suffix_doc.second, _side);
  }

 private:
  Mark *mark_;
};

template<typename AlgoCore, typename GetSuffixDocValue, typename ReportSuffixDocValue, typename IsReported>
class ComputeSuffixesByDocRMQFunctor : public ComputeSuffixesByDocFunctor {
 public:
  ComputeSuffixesByDocRMQFunctor(const AlgoCore &_core,
                                 const GetSuffixDocValue &_get_suffix_doc_value,
                                 const ReportSuffixDocValue &_report_suffix_doc_value,
                                 const IsReported &_is_reported)
      : core_{_core},
        get_suffix_doc_value_{_get_suffix_doc_value},
        report_suffix_doc_value_{_report_suffix_doc_value},
        is_reported_{_is_reported} {
  }

  template<typename AddSuffix>
  void operator()(std::size_t _sp, std::size_t _ep, AddSuffix &_add_suffix) const {
    auto report_suffix =
        [this, &_add_suffix](auto _i, const auto &_suffix_doc, OccurrenceSide _side, auto _sp, auto _ep) {
          report_suffix_doc_value_(_i, _suffix_doc, _side, _add_suffix, _sp, _ep);
        };

    GetExtremeOccurrencesRMQ<OccurrenceSide::LEFTMOST>(
        _sp, _ep, core_, core_.left_rmq, get_suffix_doc_value_, is_reported_, report_suffix);

    GetExtremeOccurrencesRMQ<OccurrenceSide::RIGHTMOST>(
        _sp, _ep, core_, core_.right_rmq, get_suffix_doc_value_, is_reported_, report_suffix);
  }

  std::vector<std::size_t> operator()(std::size_t _sp, std::size_t _ep) const override {
    std::vector<std::size_t> suffixes;
    auto add_suffix = [&suffixes](auto _suffix) {
      suffixes.emplace_back(_suffix);
    };

    (*this)(_sp, _ep, add_suffix);

    std::sort(suffixes.begin(), suffixes.end());

    return suffixes;
  }

 private:
  const AlgoCore &core_;

  const GetSuffixDocValue &get_suffix_doc_value_;
  const ReportSuffixDocValue &report_suffix_doc_value_;
  const IsReported &is_reported_;
};

template<typename AlgoCore, typename GetSuffixDocValue, typename ReportSuffixDocValue, typename IsReported>
auto MakeNewComputeSuffixesByDocRMQFunctor(const AlgoCore &_core,
                                           const GetSuffixDocValue &_get_suffix_doc_value,
                                           const ReportSuffixDocValue &_report_suffix_doc_value,
                                           const IsReported &_is_reported) {
  return new ComputeSuffixesByDocRMQFunctor<AlgoCore, GetSuffixDocValue, ReportSuffixDocValue, IsReported>(
      _core, _get_suffix_doc_value, _report_suffix_doc_value, _is_reported);
}

template<typename RAContainer, typename Report, typename DocBorderRank, typename DocBorderSelect, typename GetSuffixPosInDoc>
void ListFrequencyByDoc(const RAContainer &_suffixes,
                        Report &_report_doc_freq,
                        const DocBorderRank &_doc_border_rank,
                        const DocBorderSelect &_doc_border_select_,
                        const GetSuffixPosInDoc &_get_suffix_pos_in_doc) {
  for (int i = 0; i < _suffixes.size(); i += 2) {
    auto suffix_1 = _suffixes[i];
    auto suffix_2 = _suffixes[i + 1];
    auto doc = _doc_border_rank(suffix_1);

    if (suffix_1 == suffix_2) { // Pattern occurs exactly once
      _report_doc_freq(doc, 1);
    } else { // Pattern occurs more than once
      auto values = _get_suffix_pos_in_doc(doc, suffix_1, suffix_2);
      auto &doc_sp = values.first;
      auto &doc_ep = values.second;

      if (doc_sp > doc_ep) {
        std::swap(doc_sp, doc_ep);
      }

      _report_doc_freq(doc, doc_ep - doc_sp + 1);
    }
  }
}

template<typename LeftRMQ, typename RightRMQ, typename SetFirstRange, typename GetSuffixDocValue, typename ReportSuffixDocValue, typename IsMarkedDoc, typename MarkDoc, typename UnmarkDoc, typename GetSuffixPosInDoc, typename DocBorder, typename DocBorderRank = typename DocBorder::rank_1_type, typename DocBorderSelect = typename DocBorder::select_1_type>
class DocFreqIndexRMQ {
 public:
  DocFreqIndexRMQ(std::size_t _doc_cnt,
                  const LeftRMQ &_left_rmq,
                  const RightRMQ &_right_rmq,
                  const SetFirstRange &_set_first_range,
                  const GetSuffixDocValue &_get_suffix_doc_value,
                  const ReportSuffixDocValue &_report_suffix_doc_value,
                  const IsMarkedDoc &_is_marked_doc,
                  const MarkDoc &_mark_doc,
                  const UnmarkDoc &_unmark_doc,
                  const GetSuffixPosInDoc &_get_suffix_pos_in_doc,
                  const DocBorder &_doc_border,
                  const DocBorderRank &_doc_border_rank,
                  const DocBorderSelect &_doc_border_select)
      : doc_cnt_{_doc_cnt},
        left_rmq_{_left_rmq},
        right_rmq_{_right_rmq},
        set_first_range_{_set_first_range},
        get_suffix_doc_value_{_get_suffix_doc_value},
        report_suffix_doc_value_{_report_suffix_doc_value},
        is_marked_doc_{_is_marked_doc},
        mark_doc_{_mark_doc},
        unmark_doc_{_unmark_doc},
        get_suffix_pos_in_doc_{_get_suffix_pos_in_doc},
        doc_border_{_doc_border},
        doc_border_rank_{_doc_border_rank},
        doc_border_select_{_doc_border_select} {}

  template<typename ReportDocFreq>
  void ListWithFreq(std::size_t _sp, std::size_t _ep, ReportDocFreq _report_doc_freq) const {
    std::vector<std::size_t> suffixes = ComputeSuffixesRMQ(_sp,
                                                           _ep,
                                                           left_rmq_,
                                                           right_rmq_,
                                                           set_first_range_,
                                                           get_suffix_doc_value_,
                                                           report_suffix_doc_value_,
                                                           is_marked_doc_,
                                                           mark_doc_);

    auto report_doc_freq = [&_report_doc_freq, this](auto doc, auto freq) {
      _report_doc_freq(doc, freq);
      unmark_doc_(doc);
    };

    ListFrequencyByDoc(suffixes, report_doc_freq, doc_border_rank_, doc_border_select_, get_suffix_pos_in_doc_);
  }

  auto ListWithFreq(std::size_t _sp, std::size_t _ep) const {
    std::unordered_map<std::size_t, std::size_t> doc_freqs;

    auto report = [&doc_freqs](auto doc, auto cnt) {
      doc_freqs[doc] = cnt;
    };

    ListWithFreq(_sp, _ep, report);

    return doc_freqs;
  }

 private:
  std::size_t doc_cnt_;

  const LeftRMQ &left_rmq_;
  const RightRMQ &right_rmq_;
  const SetFirstRange &set_first_range_;

  const GetSuffixDocValue &get_suffix_doc_value_;
  const ReportSuffixDocValue &report_suffix_doc_value_;
  const IsMarkedDoc &is_marked_doc_;
  const MarkDoc &mark_doc_;
  const UnmarkDoc &unmark_doc_;

  const GetSuffixPosInDoc &get_suffix_pos_in_doc_;

  const DocBorder &doc_border_;
  const DocBorderRank &doc_border_rank_;
  const DocBorderSelect &doc_border_select_;
};

template<typename LeftRMQ, typename RightRMQ, typename TransformRange, typename GetSuffixDocValue, typename ReportSuffixDocValue, typename IsMarkedDoc, typename MarkDoc, typename UnmarkDoc, typename DocBorder, typename GetSuffixPosInDoc, typename DocBorderRank = typename DocBorder::rank_1_type, typename DocBorderSelect = typename DocBorder::select_1_type>
auto MakeDocFreqIndexRMQ(std::size_t _doc_cnt,
                         const LeftRMQ &_left_rmq,
                         const RightRMQ &_right_rmq,
                         const TransformRange &_transform_range,
                         const GetSuffixDocValue &_get_suffix_doc_value,
                         const ReportSuffixDocValue &_report_suffix_doc_value,
                         const IsMarkedDoc &_is_marked_doc,
                         const MarkDoc &_mark_doc,
                         const UnmarkDoc &_unmark_doc,
                         const GetSuffixPosInDoc &_get_suffix_pos_in_doc,
                         const DocBorder &_doc_border,
                         const DocBorderRank &_doc_border_rank,
                         const DocBorderSelect &_doc_border_select) {
  return DocFreqIndexRMQ<LeftRMQ,
                         RightRMQ,
                         TransformRange,
                         GetSuffixDocValue,
                         ReportSuffixDocValue,
                         IsMarkedDoc,
                         MarkDoc,
                         UnmarkDoc,
                         GetSuffixPosInDoc,
                         DocBorder,
                         DocBorderRank,
                         DocBorderSelect>(_doc_cnt,
                                          _left_rmq,
                                          _right_rmq,
                                          _transform_range,
                                          _get_suffix_doc_value,
                                          _report_suffix_doc_value,
                                          _is_marked_doc,
                                          _mark_doc,
                                          _unmark_doc,
                                          _get_suffix_pos_in_doc,
                                          _doc_border,
                                          _doc_border_rank,
                                          _doc_border_select);
}
}

#endif //DRET_DOC_FREQ_INDEX_RMQ_H_
