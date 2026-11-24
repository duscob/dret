//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 12/22/2024.
//

#pragma once

#include <sdsl/int_vector.hpp>

namespace sdsl {

inline bit_vector operator~(bit_vector _bv) {
  _bv.flip();

  return _bv;
}

inline bit_vector operator&(bit_vector _bv1, const bit_vector &_bv2) {
  _bv1 &= _bv2;

  return _bv1;
}

inline bit_vector operator|(bit_vector _bv1, const bit_vector &_bv2) {
  _bv1 |= _bv2;

  return _bv1;
}

}
