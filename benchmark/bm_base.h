//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/8/26.
//

#pragma once

#include <any>
#include <functional>
#include <ostream>
#include <vector>


using GenericStorage = std::map<std::string, std::any>;
using ExternalGenericStorage = std::reference_wrapper<GenericStorage>;

//~~~~~~~


template <typename Char>
std::ostream& operator<<(std::ostream& os, const std::vector<Char>& seq) {
  for (auto c : seq) {
    os << c << " ";
  }

  return os;
}

//~~~~~~~


template <typename Char>
std::vector<Char> base64_decode(const std::vector<Char>& seq, bool remove_linebreaks = false) {
  return seq;
}
