//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/20.
//

#ifndef DRET_DOC_FREQ_INDEX_SADA_H_
#define DRET_DOC_FREQ_INDEX_SADA_H_

#include <cstddef>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "doc_freq_index_rmq.h"


namespace dret {

template<OccurrenceSide _side, typename RMQ, typename GetValue, typename IsReported, typename Report>
void GetExtremeOccurrencesInRange(std::size_t _sp,
                                  std::size_t _ep,
                                  const RMQ &_rmq,
                                  const GetValue &_get_value,
                                  const IsReported &_is_reported,
                                  Report &_report) {
  using Range = std::pair<std::size_t, std::size_t>;
  std::stack<Range> stack;

  stack.emplace(_sp, _ep);
  while (!stack.empty()) {
    auto[rsp, rep] = stack.top();
    stack.pop();

    if (rsp <= rep) {
      std::size_t idx = _rmq(rsp, rep);

      auto value = _get_value(idx);

      if (!_is_reported(idx, value)) {
        _report(idx, value);

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

template<typename CSA,
    typename DocISA,
    typename RangeMinQuery,
    typename RangeMaxQuery,
    typename DocBorder,
    typename DocBorderRank = typename DocBorder::rank_1_type,
    typename DocBorderSelect = typename DocBorder::select_1_type>
class DocFreqIndexSada {
 public:
  DocFreqIndexSada(const CSA &_csa,
                   size_t _doc_cnt,
                   const DocISA &_doc_isa,
                   const DocBorder &_doc_border,
                   const DocBorderRank &_doc_border_rank,
                   const DocBorderSelect &_doc_border_select,
                   const RangeMinQuery &_range_min_query,
                   const RangeMaxQuery &_range_max_query)
      : csa_{_csa},
        doc_cnt_{_doc_cnt},
        doc_isa_{_doc_isa},
        doc_border_{_doc_border},
        doc_border_rank_{_doc_border_rank},
        doc_border_select_{_doc_border_select},
        range_min_query_{_range_min_query},
        range_max_query_{_range_max_query},
        marked_docs_{bit_vector(doc_cnt_, 0), bit_vector(doc_cnt_, 0)} {
  }

  template<typename Pattern, typename ReportDocFreq>
  void ListWithFreq(const Pattern &_pattern, ReportDocFreq _report_doc_freq) const {
    auto[sp, ep] = csa_.Search(_pattern);

    // Compute TF-IDF
    auto get_suffix_doc = [this](std::size_t idx) {
      auto suffix = csa_[idx];
      return std::make_pair(suffix, doc_border_rank_(suffix));
    };

    std::size_t left_or_right; // Set if we will mark the leftmost or rightmost occurrences. (0 == left; 1 == right)
    auto is_reported = [this, &left_or_right](auto _i, const auto &_suffix_doc) {
      return marked_docs_[left_or_right][_suffix_doc.second];
    };

    std::vector<std::size_t> suffixes;
    auto report_suffix = [this, &suffixes, &left_or_right](auto _i, const auto &_suffix_doc) {
      suffixes.emplace_back(_suffix_doc.first);
      marked_docs_[left_or_right][_suffix_doc.second] = 1;
    };

    left_or_right = 0;
    GetExtremeOccurrencesInRange<OccurrenceSide::LEFTMOST>(sp,
                                                           ep,
                                                           range_min_query_,
                                                           get_suffix_doc,
                                                           is_reported,
                                                           report_suffix);

    left_or_right = 1;
    GetExtremeOccurrencesInRange<OccurrenceSide::RIGHTMOST>(sp,
                                                            ep,
                                                            range_max_query_,
                                                            get_suffix_doc,
                                                            is_reported,
                                                            report_suffix);

    assert(suffixes.size() % 2 == 0);

    sort(suffixes.begin(), suffixes.end());

    for (int i = 0; i < suffixes.size(); i += 2) {
      auto suffix_1 = suffixes[i];
      auto suffix_2 = suffixes[i + 1];
      auto doc = doc_border_rank_(suffix_1);

      // Reset marked documents
      marked_docs_[0][doc] = 0; // Marked documents with leftmost occurrences
      marked_docs_[1][doc] = 0; // Marked documents with rightmost occurrences

      if (suffix_1 == suffix_2) { // Pattern occurs exactly once
        _report_doc_freq(doc, 1);
      } else { // Pattern occurs more than once
        std::size_t doc_begin = doc ? doc_border_select_(doc) + 1 : 0;
        std::size_t doc_sp = doc_isa_[doc][suffix_1 - doc_begin];
        std::size_t doc_ep = doc_isa_[doc][suffix_2 - doc_begin];

        // TODO Use absolute value
        if (doc_sp > doc_ep) {
          std::swap(doc_sp, doc_ep);
        }

        _report_doc_freq(doc, doc_ep - doc_sp + 1);
      }
    }
  }

  template<typename Pattern>
  auto ListWithFreq(const Pattern &_pattern) const {
    std::unordered_map<std::size_t, std::size_t> doc_freqs;

    auto report = [&doc_freqs](auto doc, auto cnt) {
      doc_freqs[doc] = cnt;
    };

    ListWithFreq(_pattern, report);

    return doc_freqs;
  }

 private:
  const CSA &csa_;
  const DocISA &doc_isa_;
  const DocBorder &doc_border_;
  const DocBorderRank &doc_border_rank_;
  const DocBorderSelect &doc_border_select_;

  std::size_t doc_cnt_;
  const RangeMinQuery &range_min_query_;
  const RangeMaxQuery &range_max_query_;

  using bit_vector = sdsl::bit_vector;
  mutable std::array<bit_vector, 2> marked_docs_;
};

template<typename CSA,
    typename DocISA,
    typename RangeMinQuery,
    typename RangeMaxQuery,
    typename DocBorder,
    typename DocBorderRank = typename DocBorder::rank_1_type,
    typename DocBorderSelect = typename DocBorder::select_1_type>
auto MakeDocFreqIndexSada(const CSA &_csa,
                          std::size_t _doc_cnt,
                          const DocISA &_doc_isa,
                          const DocBorder &_doc_border,
                          const DocBorderRank &_doc_border_rank,
                          const DocBorderSelect &_doc_border_select,
                          const RangeMinQuery &_range_min_query,
                          const RangeMaxQuery &_range_max_query) {
  return DocFreqIndexSada<CSA, DocISA, RangeMinQuery, RangeMaxQuery, DocBorder, DocBorderRank, DocBorderSelect>(
      _csa,
      _doc_cnt,
      _doc_isa,
      _doc_border,
      _doc_border_rank,
      _doc_border_select,
      _range_min_query,
      _range_max_query);
}

void ConstructPrevDocArray(const sdsl::int_vector<> &_da, std::size_t _doc_cnt, sdsl::int_vector<> &_prev_docs) {
  _prev_docs = sdsl::int_vector<>(_da.size(), 0, sdsl::bits::hi(_da.size()) + 1);
  sdsl::int_vector<> last_occ(_doc_cnt + 1, 0, sdsl::bits::hi(_da.size()) + 1);
  for (std::size_t i = 0; i < _da.size(); ++i) {
    std::size_t doc = _da[i];
    _prev_docs[i] = last_occ[doc];
    last_occ[doc] = i;
  }
}

void ConstructNextDocArray(const sdsl::int_vector<> &_da, std::size_t _doc_cnt, sdsl::int_vector<> &_next_docs) {
  _next_docs = sdsl::int_vector<>(_da.size(), 0, sdsl::bits::hi(_da.size()) + 1);
  sdsl::int_vector<> last_occ(_doc_cnt + 1, _da.size(), sdsl::bits::hi(_da.size()) + 1);
  for (std::size_t i = 0, j = _da.size() - 1; i < _da.size(); ++i, --j) {
    std::size_t doc = _da[j];
    _next_docs[j] = last_occ[doc];
    last_occ[doc] = j;
  }
}

}

#endif //DRET_DOC_FREQ_INDEX_SADA_H_
