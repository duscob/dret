//
// Created by dustin on 6/4/20.
//

#ifndef DRET_DOC_FREQ_INDEX_WT_H_
#define DRET_DOC_FREQ_INDEX_WT_H_

#include "doc_freq_index.h"

namespace dret {

template<typename CSA, typename WTOnDA>
class DocFreqIndexWT : public DocFreqIndex {
 public:
  DocFreqIndexWT(std::size_t _ndocs, const CSA &_csa, const WTOnDA &_wt_da)
      : ndocs_{_ndocs}, csa_{_csa}, wt_on_da_{_wt_da} {
  }

  std::unordered_map<std::size_t, std::size_t> Search(const std::string &_pattern) const override {
    auto values = csa_.Search(_pattern);
    auto &sp = values.first;
    auto &ep = values.second;

    std::unordered_map<std::size_t, std::size_t> doc_freqs;
    for (std::size_t i = 0; i < ndocs_; ++i) {
      auto c_r = wt_on_da_.rank(ep + 1, i);
      auto c_l = wt_on_da_.rank(sp, i);

      if (c_r != c_l)
        doc_freqs[i] = c_r - c_l;
    }

    return doc_freqs;
  }

 private:
  std::size_t ndocs_;
  const CSA &csa_;
  const WTOnDA &wt_on_da_;
};

template<typename CSA, typename WTOnDA>
auto MakeNewDocFreqIndexWT(std::size_t _ndocs, const CSA &_csa, const WTOnDA &_wt_da) {
  return new DocFreqIndexWT<CSA, WTOnDA>{_ndocs, _csa, _wt_da};
}
}

#endif //DRET_DOC_FREQ_INDEX_WT_H_
