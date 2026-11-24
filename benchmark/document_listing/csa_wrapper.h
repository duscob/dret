//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/23/20.
//

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

class CSAWrapper {
 public:
  virtual std::vector<std::size_t> Locate(const std::string& _pattern) const = 0;
  virtual std::pair<std::size_t, std::size_t> Search(const std::string& _pattern) const = 0;
};

template <typename RIndex>
class RIndexWrapper : public CSAWrapper {
 public:
  explicit RIndexWrapper(std::shared_ptr<RIndex> _r_index) : r_index_{_r_index} {}

  std::vector<std::size_t> Locate(const std::string& _pattern) const override { return r_index_->Locate(_pattern); }

  // template <typename Pattern>
  // auto Locate(const Pattern& _pattern) const {
  //   return const_cast<RIndex*>(r_index_)->Locate(const_cast<std::string&>(_pattern));
  // }

  std::pair<std::size_t, std::size_t> Search(const std::string& _pattern) const override {
    return r_index_->Count(_pattern);
  }

  // template <typename Pattern>
  // auto Search(const Pattern& _pattern) const {
  //   return const_cast<RIndex*>(r_index_)->Count(const_cast<std::string&>(_pattern));
  // }

 private:
  std::shared_ptr<RIndex> r_index_;
};
