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

template<typename VSLP, typename SLP>
class VirtualSLP : public VSLP {
 public:
  template<typename ..._Args>
  VirtualSLP(const SLP *_slp, _Args ..._args): slp_{_slp}, VSLP(_args...) {}

  auto Sigma() const {
    return slp_->Sigma();
  }

  std::pair<std::size_t, std::size_t> operator[](std::size_t i) const {
    if (i <= this->sigma_)
      return (*slp_)[i];
    else {
      auto v =  VSLP::operator[](i);
      return {v.first, v.second};
    }
//    return i <= this->sigma_ ? (*slp_)[i] : VSLP::operator[](i);
  }

  auto IsTerminal(std::size_t i) const {
    return slp_->IsTerminal(i);
  }

  auto SpanLength(std::size_t i) const {
    return i <= this->sigma_ ? slp_->SpanLength(i) : VSLP::SpanLength(i);
  }

 protected:
  const SLP *slp_;
};

template<typename Container>
class VirtualContainer : public Container {
 public:
  template<typename ..._Args>
  VirtualContainer(const Container *_container, _Args ..._args): container_{_container}, Container(_args...) {}

  auto &operator[](std::size_t i) {
    return i < container_->size() ? (*const_cast<Container *>(container_))[i] : Container::operator[](
        i - container_->size());
  }

  const auto &operator[](std::size_t i) const {
    return i < container_->size() ? (*container_)[i] : Container::operator[](i - container_->size());
  }

 protected:
  const Container *container_;
};

template<typename SLP, typename Root, typename BVs, typename ReportFirstOcc, typename ReportLastOcc>
void ReportAllFirstLastOccs(const SLP &_slp,
                            const Root &_root,
                            const BVs &_f_occs_bvs,
                            const BVs &_l_occs_bvs,
                            ReportFirstOcc &_report_f_occ,
                            ReportLastOcc &_report_l_occ) {
  std::remove_cv_t<std::remove_reference_t<decltype(_f_occs_bvs[1])>> terms(_slp.Sigma() + 1, 1);
  terms &= ~_f_occs_bvs[_root] | _l_occs_bvs[_root];

  std::size_t i = 0;
  for (const auto &term : terms) {
    if (term) {
      _report_f_occ(i, FindFirstOcc(_slp, _root, i, _f_occs_bvs, _l_occs_bvs));

      _report_l_occ(i, FindLastOcc(_slp, _root, i, _f_occs_bvs, _l_occs_bvs));
    }

    ++i;
  }
}

template<typename VSLP, typename SLP, typename Cover, typename BVs, typename ReportFirstOcc, typename ReportLastOcc>
void FindAllFirstLastOccs(const SLP &_slp,
                          const Cover &_cover,
                          const BVs &_f_occs_bvs,
                          const BVs &_l_occs_bvs,
                          ReportFirstOcc &_report_f_occ,
                          ReportLastOcc &_report_l_occ) {
  if (_cover.size() == 1) {
    ReportAllFirstLastOccs(_slp, *_cover.begin(), _f_occs_bvs, _l_occs_bvs, _report_f_occ, _report_l_occ);
    return;
  }

  VirtualSLP<VSLP, SLP> vslp(&_slp, _slp.Variables());

  auto build_vslp = [&vslp](auto id, auto left, auto right) {
    vslp.AddRule(left, right, vslp.SpanLength(left) + vslp.SpanLength(right));
  };
  BuildBinaryTree(begin(_cover), end(_cover), _slp.Variables() + 1, build_vslp);

  VirtualContainer<BVs> v_f_occs_bvs(&_f_occs_bvs, vslp.Variables() - _slp.Variables());
  VirtualContainer<BVs> v_l_occs_bvs(&_l_occs_bvs, vslp.Variables() - _slp.Variables());

  BuildFirstAndLastOccs(vslp, vslp.Start(), v_f_occs_bvs, v_l_occs_bvs);

  ReportAllFirstLastOccs(vslp, vslp.Start(), v_f_occs_bvs, v_l_occs_bvs, _report_f_occ, _report_l_occ);
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
