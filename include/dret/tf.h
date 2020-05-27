//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/30/19.
//

#ifndef DRET_TF_H_
#define DRET_TF_H_

#include <cstddef>
#include <utility>
#include <cassert>

#include "algorithm.h"

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
  } else {
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
 * Find the relative position of the first occurrence for the terminal (_term) in the expansion of the non-terminal (_var).
 *
 * The algorithm uses the bitvectors that mark the first and last occurrences of each terminal.
 *
 * @tparam _SLP Straight-line program
 * @tparam _BVs Container of bitvectors
 *
 * @param _slp
 * @param _var
 * @param _term
 * @param _f_occs_bvs bitvectors marked with the first (leftmost) occurrences
 * @param _l_occs_bvs bitvectors marked with the last (rightmost) occurrences
 *
 * @return relative position (>= 1) or 0 (if terminal doesn't appear)
 */
template<typename _SLP, typename _BVs>
std::size_t FindFirstOcc(const _SLP &_slp,
                         std::size_t _var,
                         std::size_t _term,
                         const _BVs &_f_occs_bvs,
                         const _BVs &_l_occs_bvs) {
  assert(_term <= _slp.Sigma());

  if (_slp.IsTerminal(_var)) {
    return _var == _term ? 1 : 0;
  }

  auto f_bit = _f_occs_bvs[_var][_term];
  auto l_bit = _l_occs_bvs[_var][_term];

  if (f_bit & ~l_bit) {
    return 0;
  }

  auto children = _slp[_var];
  if (!f_bit) {
    return FindFirstOcc(_slp, children.first, _term, _f_occs_bvs, _l_occs_bvs);
  } else {
    return _slp.SpanLength(children.first) + FindFirstOcc(_slp, children.second, _term, _f_occs_bvs, _l_occs_bvs);
  }
}

/**
 * Find the relative position of the last occurrence for the terminal (_term) in the expansion of the non-terminal (_var).
 *
 * The algorithm uses the bitvectors that mark the first and last occurrences of each terminal.
 *
 * @tparam _SLP Straight-line program
 * @tparam _BVs Container of bitvectors
 *
 * @param _slp
 * @param _var
 * @param _term
 * @param _f_occs_bvs bitvectors marked with the first (leftmost) occurrences
 * @param _l_occs_bvs bitvectors marked with the last (rightmost) occurrences
 *
 * @return relative position (>= 1) or 0 (if terminal doesn't appear)
 */
template<typename _SLP, typename _BVs>
std::size_t FindLastOcc(const _SLP &_slp,
                        std::size_t _var,
                        std::size_t _term,
                        const _BVs &_f_occs_bvs,
                        const _BVs &_l_occs_bvs) {
  assert(_term <= _slp.Sigma());

  if (_slp.IsTerminal(_var)) {
    return _var == _term ? 1 : 0;
  }

  auto f_bit = _f_occs_bvs[_var][_term];
  auto l_bit = _l_occs_bvs[_var][_term];

  if (f_bit & ~l_bit) {
    return 0;
  }

  auto children = _slp[_var];
  if (l_bit) {
    return _slp.SpanLength(children.first) + FindLastOcc(_slp, children.second, _term, _f_occs_bvs, _l_occs_bvs);
  } else {
    return FindLastOcc(_slp, children.first, _term, _f_occs_bvs, _l_occs_bvs);
  }
}

template<typename _SLP>
class VirtualSLP : public _SLP {
 public:
  template<typename ..._Args>
  VirtualSLP(const _SLP *_slp, _Args ..._args): slp_{_slp}, _SLP(_args...) {}

  auto Sigma() const {
    return slp_->Sigma();
  }

  auto operator[](std::size_t i) const {
    return i <= this->sigma_ ? (*slp_)[i] : _SLP::operator[](i);
  }

  auto IsTerminal(std::size_t i) const {
    return slp_->IsTerminal(i);
  }

  auto SpanLength(std::size_t i) const {
    return i <= this->sigma_ ? slp_->SpanLength(i) : _SLP::SpanLength(i);
  }

 protected:
  const _SLP *slp_;
};

template<typename _Container>
class VirtualContainer : public _Container {
 public:
  template<typename ..._Args>
  VirtualContainer(const _Container *_container, _Args ..._args): container_{_container}, _Container(_args...) {}

  auto &operator[](std::size_t i) {
    return i < container_->size() ? (*const_cast<_Container *>(container_))[i] : _Container::operator[](
        i - container_->size());
  }

  const auto &operator[](std::size_t i) const {
    return i < container_->size() ? (*container_)[i] : _Container::operator[](i - container_->size());
  }

 protected:
  const _Container *container_;
};

template<typename _SLP, typename _Cover, typename _BVs, typename _ReportFirstOcc, typename _ReportLastOcc>
void FindAllFirstLastOccs(const _SLP &_slp,
                          const _Cover &_cover,
                          const _BVs &_f_occs_bvs,
                          const _BVs &_l_occs_bvs,
                          _ReportFirstOcc &_report_f_occ,
                          _ReportLastOcc &_report_l_occ) {
  VirtualSLP<_SLP> vslp(&_slp, _slp.Variables());

  auto build_vslp = [&vslp](auto id, auto left, auto right) {
    vslp.AddRule(left, right, vslp.SpanLength(left) + vslp.SpanLength(right));
  };
  BuildBinaryTree(begin(_cover), end(_cover), _slp.Variables() + 1, build_vslp);

  VirtualContainer<_BVs> v_f_occs_bvs(&_f_occs_bvs, vslp.Variables() - _slp.Variables());
  VirtualContainer<_BVs> v_l_occs_bvs(&_l_occs_bvs, vslp.Variables() - _slp.Variables());

  BuildFirstAndLastOccs(vslp, vslp.Start(), v_f_occs_bvs, v_l_occs_bvs);

  std::remove_cv_t<std::remove_reference_t<decltype(_f_occs_bvs[1])>> terms(_slp.Sigma() + 1, 1);
  terms &= ~v_f_occs_bvs[vslp.Start()] | v_l_occs_bvs[vslp.Start()];

  std::size_t i = 0;
  for (const auto &term : terms) {
    if (term) {
      _report_f_occ(i, FindFirstOcc(vslp, vslp.Start(), i, v_f_occs_bvs, v_l_occs_bvs));

      _report_l_occ(i, FindLastOcc(vslp, vslp.Start(), i, v_f_occs_bvs, v_l_occs_bvs));
    }

    ++i;
  }
}

sdsl::bit_vector operator~(sdsl::bit_vector _bv) {
  _bv.flip();

  return _bv;
}

sdsl::bit_vector operator&(sdsl::bit_vector _bv1, const sdsl::bit_vector &_bv2) {
  _bv1 &= _bv2;

  return _bv1;
}

sdsl::bit_vector operator|(sdsl::bit_vector _bv1, const sdsl::bit_vector &_bv2) {
  _bv1 |= _bv2;

  return _bv1;
}
}
#endif //DRET_TF_H_
