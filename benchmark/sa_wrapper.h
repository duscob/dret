//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/24/20.
//

#ifndef DRET_BENCHMARK_SA_WRAPPER_H_
#define DRET_BENCHMARK_SA_WRAPPER_H_

#include <utility>
#include <vector>
#include <cstddef>
#include <grammar/differential_slp.h>

class SAWrapper {
 public:
  using DocSuffix = std::pair<std::size_t, std::size_t>; // {doc; suffix}
  typedef vector<DocSuffix> Suffixes;

  virtual DocSuffix operator()(std::size_t _i) const {
    return (*this)(_i, 1)[0];
  }

  virtual Suffixes operator()(std::size_t _i, std::size_t _n) const = 0;
};

template<typename DSLP, typename GetDoc>
class DSAWrapper : public SAWrapper {
 public:
  DSAWrapper(const DSLP *_dslp, const GetDoc *_get_doc) : dslp_{_dslp}, get_doc_{_get_doc} {}

  Suffixes operator()(std::size_t _i, std::size_t _n) const override {
    Suffixes values;
    values.reserve(_n);

    auto report = [this, &values](std::size_t _suffix) {
      values.emplace_back((*get_doc_)(_suffix), _suffix);
    };

    grammar::ExpandDifferentialSLP(*dslp_, _i, _i + _n - 1, report);

    return values;
  }

 private:
  const DSLP *dslp_;
  const GetDoc *get_doc_;
};

template<typename WTOnDA>
class WTOnDAWrapper : public SAWrapper {
 public:
  WTOnDAWrapper(const WTOnDA *_wt_da) : wt_da_{_wt_da} {}

  Suffixes operator()(std::size_t _i, std::size_t _n) const override {
    Suffixes values;
    values.reserve(_n);

    for (int i = _i; i < _i + _n; ++i) {
      auto v = wt_da_->inverse_select(i);

      values.emplace_back(v.second, v.first);
    }

    return values;
  }
 private:
  const WTOnDA *wt_da_;
};

#endif //DRET_BENCHMARK_SA_WRAPPER_H_
