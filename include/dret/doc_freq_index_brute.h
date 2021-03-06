//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/25/20.
//

#ifndef DRET_DOC_FREQ_INDEX_BRUTE_H_
#define DRET_DOC_FREQ_INDEX_BRUTE_H_

#include <cstddef>
#include <unordered_map>

#include "doc_freq_index.h"

namespace dret {

template<typename CSA, typename GetDoc, typename Pattern, typename CountDoc>
void ListFrequencyByDoc(const CSA &_csa, const GetDoc &_get_doc, const Pattern &_pattern, CountDoc _count_doc) {
  auto occurrences = _csa.Locate(_pattern);

  for (const auto &occ : occurrences) {
    auto doc = _get_doc(occ);
    _count_doc(doc);
  }
}

template<typename CSA, typename GetDoc>
class DocFreqIndexBrute : public DocFreqIndex {
 public:
  DocFreqIndexBrute(const CSA &_csa, const GetDoc &_get_doc) : csa_{_csa}, get_doc_{_get_doc} {}

  std::unordered_map<std::size_t, std::size_t> Search(const std::string &_pattern) const override {
    return ListWithFreq(_pattern);
  }

  template<typename Pattern>
  auto ListWithFreq(const Pattern &_pattern) const {
    std::unordered_map<std::size_t, std::size_t> doc_freqs;

    auto add_doc = [&doc_freqs](auto doc) {
      ++doc_freqs[doc];
    };

    ListFrequencyByDoc(csa_, get_doc_, _pattern, add_doc);

    return doc_freqs;
  }

  template<typename Pattern, typename ReportDocFreq>
  void ListWithFreq(const Pattern &_pattern, const ReportDocFreq &_report) const {
    auto doc_freqs = ListWithFreq(_pattern);

    for (const auto &item : doc_freqs) {
      _report(item.first, item.second);
    }
  }

 private:
  const CSA &csa_;
  const GetDoc &get_doc_;
};

template<typename CSA, typename GetDoc>
auto MakeDocFreqIndexBrute(const CSA &_csa, const GetDoc &_get_doc) {
  return DocFreqIndexBrute<CSA, GetDoc>{_csa, _get_doc};
}

template<typename CSA, typename GetDoc>
auto MakeNewDocFreqIndexBrute(const CSA &_csa, const GetDoc &_get_doc) {
  return new DocFreqIndexBrute<CSA, GetDoc>{_csa, _get_doc};
}

}

#endif //DRET_DOC_FREQ_INDEX_BRUTE_H_
