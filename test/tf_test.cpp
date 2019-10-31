//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/30/19.
//

#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/bit_vectors.hpp>

#include <grammar/slp.h>

#include "dret/tf.h"


using RightHand = std::pair<std::size_t, std::size_t>;
using Rules = std::vector<RightHand>;
using Occs = std::vector<sdsl::bit_vector>;
using FLOccs = std::pair<Occs, Occs>;


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


class FirstLastOccs_Test : public testing::TestWithParam<std::tuple<std::size_t, Rules, FLOccs>> {
 protected:
  grammar::SLP<> slp_;

  void SetUp() override {
    auto &sigma = std::get<0>(GetParam());
    auto &rules = std::get<1>(GetParam());

    slp_.Reset(sigma);
    for (auto &&rule : rules) {
      slp_.AddRule(rule.first, rule.second);
    }
  }
};


TEST_P(FirstLastOccs_Test, Build) {
  std::pair<Occs, Occs> occs = dret::BuildFirstAndLastOccs<sdsl::bit_vector>(slp_);

  auto &e_f_occs = std::get<2>(GetParam()).first;
  EXPECT_THAT(occs.first, testing::ElementsAreArray(e_f_occs));

  auto &e_l_occs = std::get<2>(GetParam()).second;
  EXPECT_THAT(occs.second, testing::ElementsAreArray(e_l_occs));
}


INSTANTIATE_TEST_SUITE_P(FirstLastOccs,
                         FirstLastOccs_Test,
                         testing::Values(
                             std::make_tuple(2ul,
                                             Rules{{1, 2}, {3, 3}},
                                             FLOccs{
                                                 {sdsl::bit_vector(), {1, 0, 1}, {1, 1, 0}, {1, 0, 1}, {1, 0, 0}},
                                                 {sdsl::bit_vector(), {0, 1, 0}, {0, 0, 1}, {0, 0, 1}, {0, 1, 1}}}),
                             std::make_tuple(4ul,
                                             Rules{{2, 1}, {3, 5}, {3, 3}, {2, 5}, {4, 6}, {8, 1}, {6, 7}, {11, 6},
                                                   {9, 12}, {13, 10}},
                                             FLOccs{
                                                 {sdsl::bit_vector(), {1, 0, 1, 1, 1}, {1, 1, 0, 1, 1}, {1, 1, 1, 0, 1}, {1, 1, 1, 1, 0}, {1, 1, 0, 1, 1}, {1, 1, 1, 0, 1}, {1, 1, 1, 0, 1}, {1, 1, 0, 1, 1}, {1, 1, 1, 1, 0}, {1, 0, 0, 1, 1}, {1, 0, 0, 0, 1}, {1, 0, 0, 0, 1}, {1, 0, 0, 0, 0}, {1, 0, 0, 0, 0}},
                                                 {sdsl::bit_vector(), {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 0, 0, 0}, {0, 1, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 1, 1, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 0, 0, 0}, {0, 0, 0, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 0, 0}}})/*,
                             std::make_tuple(3ul,
                                             Rules{{2, 3}, {4, 1}, {5, 3}, {2, 1}, {1, 3}, {1, 4}, {7, 2}, {1, 5},
                                                   {9, 4}, {11, 8}, {12, 5}, {6, 10}, {15, 6},
                                                   {14, 16}, {17, 13}},
                                             FLOccs{
                                                 {},
                                                 {}})*/
                         )
);

