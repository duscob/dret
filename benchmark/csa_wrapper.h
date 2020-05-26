//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/23/20.
//

#ifndef DRET_BENCHMARK_CSA_WRAPPER_H_
#define DRET_BENCHMARK_CSA_WRAPPER_H_

#include <cstddef>
#include <vector>
#include <utility>

class CSAWrapper {
 public:
  virtual std::vector<std::size_t> Locate(const std::string &_pattern) const = 0;
  virtual std::pair<std::size_t, std::size_t> Search(const std::string &_pattern) const = 0;
};

template<typename RIndex>
class RIndexWrapper: public CSAWrapper {
 public:
  explicit RIndexWrapper(const RIndex *_r_index) : r_index_{_r_index} {}

  std::vector<std::size_t> Locate(const std::string &_pattern) const override {
    return Locate<std::string>(_pattern);
  }

  template<typename Pattern>
  auto Locate(const Pattern &_pattern) const {
    return const_cast<RIndex *>(r_index_)->locate_all(const_cast<std::string &>(_pattern));
  }

  std::pair<std::size_t, std::size_t> Search(const std::string &_pattern) const override {
    return Search<std::string>(_pattern);
  }

  template<typename Pattern>
  auto Search(const Pattern &_pattern) const {
    return const_cast<RIndex *>(r_index_)->count(const_cast<std::string &>(_pattern));
  }

 private :
  const RIndex *r_index_;
};

#endif //DRET_BENCHMARK_CSA_WRAPPER_H_
