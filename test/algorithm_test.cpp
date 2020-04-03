//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/6/19.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "dret/algorithm.h"


using Seq = std::vector<std::size_t>;
using Rules = std::map<std::size_t, std::pair<std::size_t, std::size_t>>;


class BinaryTree_TF : public testing::TestWithParam<std::tuple<Seq, Rules>> {};


TEST_P(BinaryTree_TF, Build) {
  auto &seq = std::get<0>(GetParam());

  Rules rules;
  auto report = [&rules](auto _id, auto _lchild, auto _rchild) {
    rules[_id] = {_lchild, _rchild};
  };

  dret::BuildBinaryTree(seq.begin(), seq.end(), seq.size() + 1, report);

  auto &e_rules = std::get<1>(GetParam());
  EXPECT_THAT(rules, testing::ElementsAreArray(e_rules));
}


INSTANTIATE_TEST_SUITE_P(Algorithm,
                         BinaryTree_TF,
                         testing::Values(
                             std::make_tuple(Seq{1},
                                             Rules{}),
                             std::make_tuple(Seq{1, 2},
                                             Rules{{3, {1, 2}}}),
                             std::make_tuple(Seq{1, 2, 3, 4},
                                             Rules{{5, {1, 2}}, {6, {3, 4}}, {7, {5, 6}}}),
                             std::make_tuple(Seq{1, 2, 3, 4, 5},
                                             Rules{{6, {1, 2}}, {7, {3, 4}}, {8, {6, 7}}, {9, {8, 5}}})
                         )
);