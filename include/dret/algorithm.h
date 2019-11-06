//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/6/19.
//

#ifndef DRET_ALGORITHM_H_
#define DRET_ALGORITHM_H_

namespace dret {

template<typename _II, typename _Report>
auto BuildBinaryTree(_II _first, _II _last, std::size_t _initial_internal_node_id, _Report &&_report) {
  auto nnodes = distance(_first, _last);
  auto id = _initial_internal_node_id - 1;

  std::size_t nleaves = std::exp2(std::size_t(std::log2(nnodes)));

  std::vector<std::size_t> nodes(nnodes - nleaves / 2);

  auto lit_next_level = nodes.begin();
  auto create_new_nodes = [&_report, &id, &lit_next_level](auto &__fit, auto __n) {
    for (int i = 0; i < __n / 2; ++i, ++__fit, ++lit_next_level) {
      auto left = *__fit;
      auto right = *++__fit;
      _report(++id, left, right);
      *lit_next_level = id;
    }
  };

  create_new_nodes(_first, nleaves);
  std::copy(_first, _last, lit_next_level);

  while (nnodes > 1) {
    nnodes = nnodes - nleaves / 2;
    nleaves = std::exp2(std::size_t(std::log2(nnodes)));

    auto it_curr_level = lit_next_level = nodes.begin();
    create_new_nodes(it_curr_level, nleaves);
    std::copy(it_curr_level, nodes.end(), lit_next_level);
  }
}

}
#endif //DRET_ALGORITHM_H_
