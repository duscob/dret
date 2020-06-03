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

  template<typename AddDocSuffix>
  void operator()(std::size_t _sp, std::size_t _ep, AddDocSuffix &_add_doc_suffix) const {
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

    auto report_occ = [this, &_add_doc_suffix, &_sp](auto doc, auto pos) {
      _add_doc_suffix(get_suffix_(_sp + pos - 1));
    };

    FindAllFirstLastOccs<grammar::SLP<>>(slp_, cover, f_occs_bvs_, l_occs_bvs_, report_occ, report_occ);
  }

  Suffixes operator()(std::size_t _sp, std::size_t _ep) const override {
    Suffixes suffixes;
    auto add_doc_suffix = [&suffixes](auto _doc_suffix) {
      suffixes.emplace_back(_doc_suffix);
    };

    (*this)(_sp, _ep, add_doc_suffix);

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
auto MakePtrComputeSuffixesByDocGCDAFunctor(const SLP &_slp,
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

template<typename CSA, typename ComputeCover, typename GetDocs, typename GetDocFreqs>
class DocFreqIndexGCDA : public DocFreqIndex {
 public:
  DocFreqIndexGCDA(const std::shared_ptr<CSA> &_csa,
                   const std::shared_ptr<ComputeCover> &_compute_cover,
                   const std::shared_ptr<GetDocs> &_get_docs,
                   const std::shared_ptr<GetDocFreqs> &_get_doc_freqs)
      : csa_{_csa}, compute_cover_{_compute_cover}, get_docs_{_get_docs}, get_doc_freqs_{_get_doc_freqs} {
  }

  std::unordered_map<std::size_t, std::size_t> Search(const std::string &_pattern) const override {
    auto values = csa_->Search(_pattern);
    auto _sp = values.first;
    auto _ep = values.second + 1;

    auto cover = (*compute_cover_)(_sp, _ep);

    const auto &range = cover.first;
    const auto &nodes = cover.second;

    std::unordered_map<std::size_t, std::size_t> doc_freqs;
    auto add_doc = [&doc_freqs](const auto &_d) { ++doc_freqs[_d]; };

    if (nodes.empty()) {
      (*get_docs_)(_sp, _ep, add_doc);
    } else {
      (*get_docs_)(_sp, range.first, add_doc);
      (*get_docs_)(range.second, _ep, add_doc);
    }

    auto add_doc_freq = [&doc_freqs](const auto &_d, const auto &_f) {
      doc_freqs[_d] += _f;
    };

    for (const auto &node: nodes) {
      (*get_doc_freqs_)(node, add_doc_freq);
    }

    return doc_freqs;
  }

 private:
  std::shared_ptr<CSA> csa_;
  std::shared_ptr<ComputeCover> compute_cover_;
  std::shared_ptr<GetDocs> get_docs_;
  std::shared_ptr<GetDocFreqs> get_doc_freqs_;
};

template<typename CSA, typename ComputeCover, typename GetDocs, typename GetDocFreqs>
auto MakePtrDocFreqIndexGCDA(const std::shared_ptr<CSA> &_csa,
                             const std::shared_ptr<ComputeCover> &_compute_cover,
                             const std::shared_ptr<GetDocs> &_get_docs,
                             const std::shared_ptr<GetDocFreqs> &_get_doc_freqs) {
  return new DocFreqIndexGCDA<CSA, ComputeCover, GetDocs, GetDocFreqs>(_csa, _compute_cover, _get_docs, _get_doc_freqs);
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromLeft(VarType _var, std::size_t &_length, const SLP &_slp, Report &_report) {
  if (_slp.IsTerminal(_var)) {
    _report(_var);
    --_length;
    return;
  }

  const auto &children = _slp[_var];
  ExpandSLPFromLeft(children.first, _length, _slp, _report);

  if (_length) {
    ExpandSLPFromLeft(children.second, _length, _slp, _report);
  }
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromFront(VarType _var, std::size_t _length, const SLP &_slp, Report &_report) {
  do {
    const auto &cover = _slp.Cover(_var);

    auto size = cover.size();

    for (auto it = cover.begin(); _length && size; --size, ++it) {
      ExpandSLPFromLeft(*it, _length, _slp, _report);
    }

    ++_var;
  } while (_length);
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromLeft(VarType _var, std::size_t &_skip, std::size_t &_length, const SLP &_slp, Report &_report) {
  if (_slp.IsTerminal(_var)) {
    --_skip;
    return;
  }

  const auto &children = _slp[_var];
  ExpandSLPFromLeft(children.first, _skip, _length, _slp, _report);

  if (_skip) {
    ExpandSLPFromLeft(children.second, _skip, _length, _slp, _report);
  } else if (_length) {
    ExpandSLPFromLeft(children.second, _length, _slp, _report);
  }
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromFront(VarType _var, std::size_t _skip, std::size_t _length, const SLP &_slp, Report &_report) {
  const auto &cover = _slp.Cover(_var);

  auto size = cover.size();

  for (auto it = cover.begin(); _length && size; --size, ++it) {
    if (_skip) {
      ExpandSLPFromLeft(*it, _skip, _length, _slp, _report);
    } else {
      ExpandSLPFromLeft(*it, _length, _slp, _report);
    }
  }
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromRight(VarType _var, std::size_t &_length, const SLP &_slp, Report &_report) {
  if (_slp.IsTerminal(_var)) {
    _report(_var);
    --_length;
    return;
  }

  const auto &children = _slp[_var];
  ExpandSLPFromRight(children.second, _length, _slp, _report);

  if (_length) {
    ExpandSLPFromRight(children.first, _length, _slp, _report);
  }
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromBack(VarType _var, std::size_t _length, const SLP &_slp, Report &_report) {
  const auto &cover = _slp.Cover(_var);

  auto size = cover.size();

  for (auto it = cover.end() - 1; _length && size; --size, --it) {
    ExpandSLPFromRight(*it, _length, _slp, _report);
  }
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromRight(VarType _var, std::size_t &_skip, std::size_t &_length, const SLP &_slp, Report &_report) {
  if (_slp.IsTerminal(_var)) {
    --_skip;
    return;
  }

  const auto &children = _slp[_var];
  ExpandSLPFromRight(children.second, _skip, _length, _slp, _report);

  if (_skip) {
    ExpandSLPFromRight(children.first, _skip, _length, _slp, _report);
  } else if (_length) {
    ExpandSLPFromRight(children.first, _length, _slp, _report);
  }
}

template<typename VarType, typename SLP, typename Report>
void ExpandSLPFromBack(VarType _var, std::size_t _skip, std::size_t _length, const SLP &_slp, Report &_report) {
  const auto &cover = _slp.Cover(_var);

  auto size = cover.size();

  for (auto it = cover.end() - 1; _length && size; --size, --it) {
    if (_skip) {
      ExpandSLPFromRight(*it, _skip, _length, _slp, _report);
    } else {
      ExpandSLPFromRight(*it, _length, _slp, _report);
    }
  }
}

template<typename SLP, typename Report>
void ExpandSLP(const SLP &_slp, std::size_t _bp, std::size_t _ep, Report &_report) {
  if (_bp >= _ep)
    return;

  auto leaf = _slp.Leaf(_bp);
  auto pos = _slp.Position(leaf);

  if (pos == _bp) {
    ExpandSLPFromFront(leaf, _ep - pos, _slp, _report);
  } else {
    auto next_pos = _slp.Position(leaf + 1);

    if (next_pos <= _ep) {
      ExpandSLPFromBack(leaf, next_pos - _bp, _slp, _report);

      if (next_pos < _ep) {
        ExpandSLPFromFront(leaf + 1, _ep - next_pos, _slp, _report);
      }
    } else {
      auto skip_front = _bp - pos;
      auto skip_back = next_pos - _ep;
      if (skip_front < skip_back) {
        ExpandSLPFromFront(leaf, skip_front, _ep - _bp, _slp, _report);
      } else {
        ExpandSLPFromBack(leaf, skip_back, _ep - _bp, _slp, _report);
      }
    }
  }
}

template<typename SLP>
class ExpandSLPFunctor {
 public:
  explicit ExpandSLPFunctor(const SLP &_slp) : slp_{_slp} {}

  template<typename Report>
  void operator()(std::size_t _bp, std::size_t _ep, Report &_report) const {
    ExpandSLP(slp_, _bp, _ep, _report);
  }

 private:
  const SLP &slp_;
};

template<typename SLP>
auto MakePtrExpandSLPFunctor(const SLP &_slp) {
  return new ExpandSLPFunctor<SLP>(_slp);
}

}

#endif //DRET_DOC_FREQ_INDEX_GCDA_H_
