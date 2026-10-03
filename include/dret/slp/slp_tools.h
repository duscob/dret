//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/15/26.
//

#pragma once

#include <cstddef>
#include <vector>

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

//~~~~~~~  In-order expansion with early stop  ~~~~~~~
//
// ExpandSLP above is not in position order: a range that starts in the back half
// of its sampled leaf is expanded from the leaf's end, right to left. The RMQ
// cores read a run in one pass and need its head first -- they stop there when
// the head document is already reported, and CILCP also stops after a peek at the
// second position -- so ExpandSLPUntil reports [bp, ep) strictly left to right
// and stops as soon as the report returns false.
//
// Locating bp costs what a single-position ExpandSLP(bp, bp + 1) costs: the leaf
// is entered from the end nearer to bp. From the back, the leaf's tail [bp, end)
// is collected into a buffer and replayed in order, so that case cannot stop
// early inside the leaf; it costs about the same as that single lookup.

namespace internal {

// Forward expansion of _var: skips _skip terminals, then reports up to _length.
// Returns false once the report has asked to stop.
template<typename VarType, typename SLP, typename Report>
bool ExpandUntilFromLeft(VarType _var, std::size_t &_skip, std::size_t &_length, const SLP &_slp, Report &_report) {
  if (_slp.IsTerminal(_var)) {
    if (_skip) {
      --_skip;
      return true;
    }
    --_length;
    return _report(_var);
  }

  const auto &children = _slp[_var];
  if (!ExpandUntilFromLeft(children.first, _skip, _length, _slp, _report))
    return false;
  if (_length)
    return ExpandUntilFromLeft(children.second, _skip, _length, _slp, _report);
  return true;
}

// Forward expansion of the cover of _leaf; same contract as ExpandUntilFromLeft.
template<typename VarType, typename SLP, typename Report>
bool ExpandUntilFromFront(VarType _leaf, std::size_t &_skip, std::size_t &_length, const SLP &_slp, Report &_report) {
  const auto &cover = _slp.Cover(_leaf);
  auto size = cover.size();
  for (auto it = cover.begin(); _length && size; --size, ++it) {
    if (!ExpandUntilFromLeft(*it, _skip, _length, _slp, _report))
      return false;
  }
  return true;
}

}  // namespace internal

template<typename SLP, typename Report>
void ExpandSLPUntil(const SLP &_slp, std::size_t _bp, std::size_t _ep, Report &_report) {
  if (_bp >= _ep)
    return;

  auto leaf = _slp.Leaf(_bp);
  const auto pos = _slp.Position(leaf);
  if (pos != _bp) {
    const auto next_pos = _slp.Position(leaf + 1);
    const auto end = next_pos < _ep ? next_pos : _ep;
    // Enter the leaf from whichever end is nearer to bp itself, as the
    // single-position lookup ExpandSLP(bp, bp + 1) would: the head alone decides,
    // not the length of the run behind it, since the run may stop at the head.
    const auto skip_front = _bp - pos;
    const auto skip_back = next_pos - end;
    if (skip_front < next_pos - _bp - 1) {
      std::size_t skip = skip_front, length = end - _bp;
      if (!internal::ExpandUntilFromFront(leaf, skip, length, _slp, _report))
        return;
    } else {
      thread_local std::vector<std::size_t> tail;
      tail.clear();
      auto collect = [](auto _v) { tail.emplace_back(static_cast<std::size_t>(_v)); };
      ExpandSLPFromBack(leaf, skip_back, end - _bp, _slp, collect);
      for (auto it = tail.rbegin(); it != tail.rend(); ++it) {
        if (!_report(*it))
          return;
      }
    }
    if (end == _ep)
      return;
    _bp = end;
    ++leaf;
  }

  std::size_t length = _ep - _bp;
  while (length) {
    std::size_t skip = 0;
    if (!internal::ExpandUntilFromFront(leaf, skip, length, _slp, _report))
      return;
    ++leaf;
  }
}

}