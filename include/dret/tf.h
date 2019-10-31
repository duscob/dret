//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/30/19.
//

#ifndef DRET_TF_H_
#define DRET_TF_H_

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

}
#endif //DRET_TF_H_
