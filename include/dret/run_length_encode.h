//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/31/20.
//

#ifndef DRET_RUN_LENGTH_ENCODE_H_
#define DRET_RUN_LENGTH_ENCODE_H_

#include <cstddef>
#include <iostream>
#include <iterator>

namespace dret {
/***
 * Build the Run-Length Encoding bitvector for a given sequence
 * @tparam IIter Input Iterator
 * @tparam BVRuns BitVector to store the runs of equal values
 *
 * @param _first
 * @param _last
 * @param _bv_runs
 */
template<typename IIter, typename BVRuns>
void BuildRLEncoding(IIter _first, IIter _last, BVRuns &_bv_runs) {
  if (_first == _last) return;

  _bv_runs = BVRuns(distance(_first, _last));
  std::size_t i = 0;
//  _bv_runs[i] = 0;
  _bv_runs[i] = 1;
  for (auto it = _first++; _first != _last; ++it, ++_first) {
//    _bv_runs[++i] = (*it == *_first);
    _bv_runs[++i] = (*it != *_first);
  }
}


/**
 * Run-Length Encoding
 * @tparam BVRuns BitVector to store the runs of equal values
 * @tparam BVRunsRank Rank-support over the bitvector with the runs
 */
template<typename BVRuns, typename BVRunsRank = typename BVRuns::rank_1_type>
class RLEncoding {
 public:
  RLEncoding() = default;

  template<typename _BVRuns>
  RLEncoding(const _BVRuns &_bv_runs): bv_runs_(_bv_runs), bv_runs_rank_(&bv_runs_) {}

  template<typename _BVRuns, typename _BVRunsRank>
  RLEncoding(RLEncoding<_BVRuns, _BVRunsRank> _rle): bv_runs_(_rle.GetRuns()), bv_runs_rank_(&bv_runs_) {}

  /**
   * Compute the new position on the run-length encoding sequence
   *
   * @param _pos Original Position
   * @return Run-Length Encoding Position
   */
  std::size_t operator()(std::size_t _pos) const {
//    return _pos - bv_runs_rank_(_pos + 1);
    return bv_runs_rank_(_pos + 1);
  }

  const BVRuns & GetRuns() const {
    return bv_runs_;
  }

  typedef std::size_t size_type;

  template<typename ...Args>
  size_type serialize(std::ostream &out, Args ..._args/*sdsl::structure_tree_node *v = nullptr, const std::string &name = ""*/) const {
    std::size_t written_bytes = 0;
    written_bytes += bv_runs_.serialize(out);
    written_bytes += bv_runs_rank_.serialize(out);

    return written_bytes;
  }

  void load(std::istream &in) {
    bv_runs_.load(in);
    bv_runs_rank_.load(in, &bv_runs_);
  }

 protected:
  BVRuns bv_runs_;
  BVRunsRank bv_runs_rank_;
};


template<typename Iter>
class RLEIterator : public std::iterator<std::input_iterator_tag,
                                         typename Iter::value_type,
                                         typename Iter::difference_type,
                                         typename Iter::pointer,
                                         typename Iter::reference> {
 public:
  RLEIterator(const Iter &_iter, Iter _eiter) : iter_(_iter), eiter_(_eiter) {}

  RLEIterator &operator++() {
    for (auto it = iter_++; iter_ != eiter_ && *it == *iter_; ++it, ++iter_);

    return *this;
  }

  RLEIterator operator++(int) {
    Iter first = iter_;
    for (auto it = iter_++; iter_ != eiter_ && *it == *iter_; ++it, ++iter_);

    return RLEIterator<Iter>(first, eiter_);
  }

  bool operator==(RLEIterator _other) const { return iter_ == _other.iter_; }

  bool operator!=(RLEIterator _other) const { return !(*this == _other); }

  typename Iter::reference operator*() const { return *iter_; }

  typename Iter::pointer operator->() const { return iter_.operator->(); }

 protected:
  Iter iter_;
  Iter eiter_;
};

template <typename Iter>
auto MakeRLEIterator(const Iter &_iter, const Iter _eiter) {
  return RLEIterator<Iter>(_iter, _eiter);
}

}

#endif //DRET_RUN_LENGTH_ENCODE_H_
