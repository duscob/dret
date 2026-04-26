//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/15/26.
//

#pragma once

namespace dret {

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
  ExpandSLPFunctor() = default;

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