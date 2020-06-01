//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/27/20.
//

#ifndef DRET_DOC_FREQ_INDEX_GCDA_H_
#define DRET_DOC_FREQ_INDEX_GCDA_H_

#include <cstddef>
#include <vector>

#include "doc_freq_index.h"
#include "tf.h"

namespace dret {

template<typename SLP, typename Report>
void ComputeSpanCoverFrom(const SLP &_slp, std::size_t _var, std::size_t _sp, const Report &_report) {
  if (_sp == 0) {
    _report(_var);
    return;
  }

  const auto &children = _slp[_var];
  auto left_length = _slp.SpanLength(children.first);

  if (left_length <= _sp) {
    ComputeSpanCoverFrom(_slp, children.second, _sp - left_length, _report);
    return;
  }

  ComputeSpanCoverFrom(_slp, children.first, _sp, _report);
  _report(children.second);
}

template<typename SLP, typename Report>
void ComputeSpanCoverTo(const SLP &_slp, std::size_t _var, std::size_t _ep, const Report &_report) {
  if (_slp.SpanLength(_var) <= _ep + 1) {
    _report(_var);
    return;
  }

  const auto &children = _slp[_var];
  auto left_length = _slp.SpanLength(children.first);

  if (_ep < left_length) {
    ComputeSpanCoverTo(_slp, children.first, _ep, _report);
    return;
  }

  _report(children.first);
  ComputeSpanCoverTo(_slp, children.second, _ep - left_length, _report);
}

template<typename SLP, typename Report>
void ComputeSpanCoverIn(const SLP &_slp, std::size_t _var, std::size_t _sp, std::size_t _ep, Report _report) {
  if (_sp > _ep)
    return;

  if (0 == _sp && _slp.SpanLength(_var) <= _ep + 1) {
    _report(_var);
    return;
  }

  const auto &children = _slp[_var];
  auto left_length = _slp.SpanLength(children.first);

  if (_ep < left_length) {
    ComputeSpanCoverIn(_slp, children.first, _sp, _ep, _report);
    return;
  }

  if (left_length <= _sp) {
    ComputeSpanCoverIn(_slp, children.second, _sp - left_length, _ep - left_length, _report);
    return;
  }

  ComputeSpanCoverFrom(_slp, children.first, _sp, _report);
  ComputeSpanCoverTo(_slp, children.second, _ep - left_length, _report);
}

template<typename SLP, typename Roots, typename RootHeadBV, typename OccsBVs, typename GetSuffix, typename RootHeadBVRank = typename RootHeadBV::rank_1_type, typename RootHeadBVSelect = typename RootHeadBV::select_1_type>
class ComputeSuffixesByDocGCDAFunctor : public ComputeSuffixesByDocFunctor {
 public:
  ComputeSuffixesByDocGCDAFunctor(const SLP &_slp,
                                  const Roots &_roots,
                                  const RootHeadBV &_root_head_bv,
                                  const RootHeadBVRank &_root_head_bv_rank,
                                  const RootHeadBVSelect &_root_head_bv_select,
                                  const OccsBVs &_f_occs_bvs,
                                  const OccsBVs &_l_occs_bvs,
                                  const GetSuffix &_get_suffix)
      : slp_{_slp},
        roots_{_roots},
        root_head_bv_{_root_head_bv},
        root_head_bv_rank_{_root_head_bv_rank},
        root_head_bv_select_{_root_head_bv_select},
        f_occs_bvs_{_f_occs_bvs},
        l_occs_bvs_{_l_occs_bvs},
        get_suffix_{_get_suffix} {
  }

  template<typename AddSuffix>
  void operator()(std::size_t _sp, std::size_t _ep, AddSuffix &_add_suffix) const {
    if (_ep < _sp)
      return;

    std::vector<std::size_t> cover;
    auto add_var = [&cover](const auto &var) {
      cover.emplace_back(var);
    };

    auto l = root_head_bv_rank_(_sp + 1) - 1;
    auto r = root_head_bv_rank_(_ep + 1) - 1;

    auto root_cover_sp = root_head_bv_select_(l + 1);

    if (l < r) {
      if (root_cover_sp < _sp) {
        // Prefix sequence not cover by a complete root
        ComputeSpanCoverFrom(slp_, roots_[l], _sp - root_cover_sp, add_var);
        ++l;
      }

      // Sequence cover by complete roots
      for (auto i = l; i < r; ++i) {
        cover.emplace_back(roots_[i]);
      }

      auto root_cover_ep = root_head_bv_select_(r + 1);
      if (root_cover_ep <= _ep) {
        // Suffix sequence not cover by a complete root
        ComputeSpanCoverTo(slp_, roots_[r], _ep - root_cover_ep, add_var);
      } else {
        cover.emplace_back(roots_[r]);
      }
    } else {
      ComputeSpanCoverIn(slp_, roots_[l], _sp - root_cover_sp, _ep - root_cover_sp, add_var);
    }

    auto report_occ = [this, &_add_suffix, &_sp](auto doc, auto pos) {
      _add_suffix(get_suffix_(_sp + pos - 1).first);
    };

    FindAllFirstLastOccs<grammar::SLP<>>(slp_, cover, f_occs_bvs_, l_occs_bvs_, report_occ, report_occ);
  }

  std::vector<std::size_t> operator()(std::size_t _sp, std::size_t _ep) const override {
    std::vector<std::size_t> suffixes;
    auto add_suffix = [&suffixes](auto _suffix) {
      suffixes.emplace_back(_suffix);
    };

    (*this)(_sp, _ep, add_suffix);

    std::sort(suffixes.begin(), suffixes.end());

    return suffixes;
  }

 private:
  const SLP &slp_;
  const Roots &roots_;
  const RootHeadBV &root_head_bv_;
  const RootHeadBVRank &root_head_bv_rank_;
  const RootHeadBVSelect &root_head_bv_select_;

  const OccsBVs &f_occs_bvs_;
  const OccsBVs &l_occs_bvs_;

  const GetSuffix &get_suffix_;
};

template<typename SLP, typename Roots, typename RootHeadBV, typename OccsBVs, typename GetSuffix, typename RootHeadBVRank = typename RootHeadBV::rank_1_type, typename RootHeadBVSelect = typename RootHeadBV::select_1_type>
auto MakeNewComputeSuffixesByDocGCDAFunctor(const SLP &_slp,
                                            const Roots &_roots,
                                            const RootHeadBV &_root_head_bv,
                                            const RootHeadBVRank &_root_head_bv_rank,
                                            const RootHeadBVSelect &_root_head_bv_select,
                                            const OccsBVs &_f_occs_bvs,
                                            const OccsBVs &_l_occs_bvs,
                                            const GetSuffix &_get_suffix) {
  return new ComputeSuffixesByDocGCDAFunctor<SLP,
                                             Roots,
                                             RootHeadBV,
                                             OccsBVs,
                                             GetSuffix,
                                             RootHeadBVRank,
                                             RootHeadBVSelect>
      (_slp, _roots, _root_head_bv, _root_head_bv_rank, _root_head_bv_select, _f_occs_bvs, _l_occs_bvs, _get_suffix);
}

}

#endif //DRET_DOC_FREQ_INDEX_GCDA_H_
