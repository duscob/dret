//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/25/20.
//

#ifndef DRET_DF_BRUTE_H_
#define DRET_DF_BRUTE_H_

#include <cstddef>
#include <unordered_map>

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
class DFIdxBrute {
 public:
  DFIdxBrute(const CSA &_csa, const GetDoc &_get_doc) : csa_{_csa}, get_doc_{_get_doc} {
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
auto MakeDFIdxBrute(const CSA &_csa, const GetDoc &_get_doc) {
  return DFIdxBrute<CSA, GetDoc>{_csa, _get_doc};
}

}

#endif //DRET_DF_BRUTE_H_
