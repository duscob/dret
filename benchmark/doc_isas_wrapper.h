//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/24/20.
//

#ifndef DRET_BENCHMARK_DOC_ISAS_WRAPPER_H_
#define DRET_BENCHMARK_DOC_ISAS_WRAPPER_H_

#include <cstddef>
#include <utility>
#include <vector>

class DocISAsWrapper {
 public:
  virtual std::pair<std::size_t, std::size_t> operator()(std::size_t _doc,
                                                         std::size_t _suffix_1,
                                                         std::size_t suffix_2) const = 0;
};

template<typename DSLP>
class DocDISAsWrapper: public DocISAsWrapper {
 public:
  DocDISAsWrapper(const DSLP *_dslp) : dslp_{_dslp} {}

  std::pair<std::size_t, std::size_t> operator()(std::size_t _doc,
                                                 std::size_t _suffix_1,
                                                 std::size_t _suffix_2) const override {
    std::vector<std::size_t> values;
    values.reserve(2);

    auto report = [&values](std::size_t _value) {
      values.emplace_back(_value);
    };

    grammar::ExpandDifferentialSLP(*dslp_, _suffix_1, _suffix_1, report);
    grammar::ExpandDifferentialSLP(*dslp_, _suffix_2, _suffix_2, report);

    return std::make_pair(values[0], values[1]);
  }

 private:
  const DSLP *dslp_;
};

#endif //DRET_BENCHMARK_DOC_ISAS_WRAPPER_H_
