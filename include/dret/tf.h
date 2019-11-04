//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/30/19.
//

#ifndef DRET_TF_H_
#define DRET_TF_H_

#include <cstddef>
#include <utility>
#include <cassert>


namespace dret {

template<typename _SLP, typename _OccBvs>
void BuildFirstAndLastOccs(const _SLP &_slp, std::size_t _var, _OccBvs &_f_occs, _OccBvs &_l_occs) {
  if (!_l_occs[_var].empty() && !_f_occs[_var].empty())
    return;

  auto n_terminals = _slp.Sigma() + 1;
  if (_slp.IsTerminal(_var)) {
    _f_occs[_var] = typename _OccBvs::value_type(n_terminals, 1);
    _f_occs[_var][_var] = 0;

    _l_occs[_var] = typename _OccBvs::value_type(n_terminals, 0);
    _l_occs[_var][_var] = 1;
  }
  else {
    auto children = _slp[_var];

    BuildFirstAndLastOccs(_slp, children.first, _f_occs, _l_occs);
    BuildFirstAndLastOccs(_slp, children.second, _f_occs, _l_occs);

    _f_occs[_var] = _f_occs[children.first] & ~_l_occs[children.first];

    _l_occs[_var] = ~_f_occs[children.second] | _l_occs[children.second];
  }
}


/**
 * Build bitvectors marking the first(leftmost) and last(rightmost) occurrences of each terminal for each non-terminal.
 *
 * @tparam _Bitvector
 * @tparam _SLP
 * @param _slp
 * @return
 */
template<typename _Bitvector, typename _SLP>
auto BuildFirstAndLastOccs(const _SLP &_slp) {
  using OccsBV = std::vector<_Bitvector>;
  auto occs = std::make_pair(OccsBV(_slp.Variables() + 1), OccsBV(_slp.Variables() + 1));

  BuildFirstAndLastOccs(_slp, _slp.Start(), occs.first, occs.second);

  return occs;
}


/**
 * Find the relative position of the first occurrence for the document (_doc) in the expansion of the non-terminal (_var).
 *
 * The algorithm uses the bitvectors that mark the first and last occurrences of each document.
 *
 * @tparam _SLP Straight-line program
 * @tparam _BVs Container of bitvectors
 *
 * @param _slp
 * @param _var
 * @param _doc
 * @param _f_occs_bvs bitvectors marked with the first (leftmost) occurrences
 * @param _l_occs_bvs bitvectors marked with the last (rightmost) occurrences
 *
 * @return relative position (>= 1) or 0 (if document doesn't appear)
 */
template<typename _SLP, typename _BVs>
std::size_t FindFirstOcc(const _SLP &_slp, std::size_t _var, std::size_t _doc, const _BVs &_f_occs_bvs, const _BVs &_l_occs_bvs) {
  assert(_doc <= _slp.Sigma());

  if (_slp.IsTerminal(_var) && _var == _doc) {
    return 1;
  }

  auto f_bit = _f_occs_bvs[_var][_doc];
  auto l_bit = _l_occs_bvs[_var][_doc];

  if (f_bit & ~l_bit) {
    return 0;
  }

  auto children = _slp[_var];
  if (!f_bit) {
    return FindFirstOcc(_slp, children.first, _doc, _f_occs_bvs, _l_occs_bvs);
  }
  else {
    return _slp.SpanLength(children.first) + FindFirstOcc(_slp, children.second, _doc, _f_occs_bvs, _l_occs_bvs);
  }
}


/**
 * Find the relative position of the last occurrence for the document (_doc) in the expansion of the non-terminal (_var).
 *
 * The algorithm uses the bitvectors that mark the first and last occurrences of each document.
 *
 * @tparam _SLP Straight-line program
 * @tparam _BVs Container of bitvectors
 *
 * @param _slp
 * @param _var
 * @param _doc
 * @param _f_occs_bvs bitvectors marked with the first (leftmost) occurrences
 * @param _l_occs_bvs bitvectors marked with the last (rightmost) occurrences
 *
 * @return relative position (>= 1) or 0 (if document doesn't appear)
 */
template<typename _SLP, typename _BVs>
std::size_t FindLastOcc(const _SLP &_slp, std::size_t _var, std::size_t _doc, const _BVs &_f_occs_bvs, const _BVs &_l_occs_bvs) {
  assert(_doc <= _slp.Sigma());

  if (_slp.IsTerminal(_var) && _var == _doc) {
    return 1;
  }

  auto f_bit = _f_occs_bvs[_var][_doc];
  auto l_bit = _l_occs_bvs[_var][_doc];

  if (f_bit & ~l_bit) {
    return 0;
  }

  auto children = _slp[_var];
  if (l_bit) {
    return _slp.SpanLength(children.first) + FindLastOcc(_slp, children.second, _doc, _f_occs_bvs, _l_occs_bvs);
  }
  else {
    return FindLastOcc(_slp, children.first, _doc, _f_occs_bvs, _l_occs_bvs);
  }
}
}
#endif //DRET_TF_H_
