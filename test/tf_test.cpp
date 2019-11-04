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
using OccsBVs = std::vector<sdsl::bit_vector>;
using FLOccsBVs = std::pair<OccsBVs, OccsBVs>;
using Occ = std::pair<std::pair<std::size_t, std::vector<std::size_t>>,
                      std::vector<std::size_t>>; // (non-terminal, docs) = relative_positions
using Occs = std::vector<Occ>;


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


template<typename ... _Args>
class TestWithSLP : public testing::TestWithParam<std::tuple<std::size_t, Rules, _Args...>> {
 protected:
  grammar::SLP<> slp_;

  void SetUp() override {
    auto &sigma = std::get<0>(this->GetParam());
    auto &rules = std::get<1>(this->GetParam());

    slp_.Reset(sigma);
    for (auto &&rule : rules) {
      slp_.AddRule(rule.first, rule.second);
    }
  }
};


class FirstLastOccs_Test : public TestWithSLP<FLOccsBVs> {
};


TEST_P(FirstLastOccs_Test, Build) {
  std::pair<OccsBVs, OccsBVs> occs_bvs = dret::BuildFirstAndLastOccs<sdsl::bit_vector>(slp_);

  auto &e_f_occs_bv = std::get<2>(GetParam()).first;
  EXPECT_THAT(occs_bvs.first, testing::ElementsAreArray(e_f_occs_bv));

  auto &e_l_occs_bv = std::get<2>(GetParam()).second;
  EXPECT_THAT(occs_bvs.second, testing::ElementsAreArray(e_l_occs_bv));
}


INSTANTIATE_TEST_SUITE_P(FirstLastOccs,
                         FirstLastOccs_Test,
                         testing::Values(
                             std::make_tuple(2ul,
                                             Rules{{1, 2}, {3, 3}},
                                             FLOccsBVs{
                                                 {sdsl::bit_vector(), {1, 0, 1}, {1, 1, 0}, {1, 0, 1}, {1, 0, 0}},
                                                 {sdsl::bit_vector(), {0, 1, 0}, {0, 0, 1}, {0, 0, 1}, {0, 1, 1}}}),
                             std::make_tuple(4ul,
                                             Rules{{2, 1}, {3, 5}, {3, 3}, {2, 5}, {4, 6}, {8, 1}, {6, 7}, {11, 6},
                                                   {9, 12}, {13, 10}},
                                             FLOccsBVs{
                                                 {sdsl::bit_vector(), {1, 0, 1, 1, 1}, {1, 1, 0, 1, 1}, {1, 1, 1, 0, 1},
                                                  {1, 1, 1, 1, 0}, {1, 1, 0, 1, 1}, {1, 1, 1, 0, 1}, {1, 1, 1, 0, 1},
                                                  {1, 1, 0, 1, 1}, {1, 1, 1, 1, 0}, {1, 0, 0, 1, 1}, {1, 0, 0, 0, 1},
                                                  {1, 0, 0, 0, 1}, {1, 0, 0, 0, 0}, {1, 0, 0, 0, 0}},
                                                 {sdsl::bit_vector(), {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0},
                                                  {0, 0, 0, 0, 1}, {0, 1, 0, 0, 0}, {0, 1, 1, 0, 0}, {0, 0, 0, 1, 0},
                                                  {0, 1, 1, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 0, 0, 0}, {0, 0, 0, 1, 0},
                                                  {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 0, 0}}})/*,
                             std::make_tuple(3ul,
                                             Rules{{2, 3}, {4, 1}, {5, 3}, {2, 1}, {1, 3}, {1, 4}, {7, 2}, {1, 5},
                                                   {9, 4}, {11, 8}, {12, 5}, {6, 10}, {15, 6},
                                                   {14, 16}, {17, 13}},
                                             FLOccsBVs{
                                                 {},
                                                 {}})*/
                         )
);


class FirstOccs_Test : public TestWithSLP<Occs> {};


TEST_P(FirstOccs_Test, Find) {
  auto[f_occs_bvs, l_occs_bvs] = dret::BuildFirstAndLastOccs<sdsl::bit_vector>(slp_);

  auto occs = std::get<2>(GetParam());

  for (const auto &occ : occs) {
    auto &nonterminal = occ.first.first;
    auto &docs = occ.first.second;
    auto &e_poss = occ.second;

    for (std::size_t i = 0; i < docs.size(); ++i) {
      auto pos = dret::FindFirstOcc(slp_, nonterminal, docs[i], f_occs_bvs, l_occs_bvs);
      EXPECT_EQ(pos, e_poss[i]) << nonterminal << ":" << docs[i];
    }
  }
}


INSTANTIATE_TEST_SUITE_P(FirstOccs,
                         FirstOccs_Test,
                         testing::Values(
                             std::make_tuple(2ul,
                                             Rules{{1, 2}, {3, 3}},
                                             Occs{{{1, {1, 2}}, {1, 0}}, {{2, {2, 1}}, {1, 0}}, {{3, {1, 2}}, {1, 2}},
                                                  {{4, {1, 2}}, {1, 2}}}),
                             std::make_tuple(4ul,
                                             Rules{{2, 1}, {3, 5}, {3, 3}, {2, 5}, {4, 6}, {8, 1}, {6, 7}, {11, 6},
                                                   {9, 12}, {13, 10}},
                                             Occs{{{1, {1, 2, 3, 4}}, {1, 0, 0, 0}}, {{2, {1, 2, 3, 4}}, {0, 1, 0, 0}},
                                                  {{3, {1, 2, 3, 4}}, {0, 0, 1, 0}}, {{4, {1, 2, 3, 4}}, {0, 0, 0, 1}},
                                                  {{5, {1, 2}}, {2, 1}}, {{6, {1, 2, 3}}, {3, 2, 1}},
                                                  {{7, {1, 3}}, {0, 1}}, {{10, {1, 2}}, {3, 1}}})
                         )
);


class LastOccs_Test : public TestWithSLP<Occs> {};


TEST_P(LastOccs_Test, Find) {
  auto[f_occs_bvs, l_occs_bvs] = dret::BuildFirstAndLastOccs<sdsl::bit_vector>(slp_);

  auto occs = std::get<2>(GetParam());

  for (const auto &occ : occs) {
    auto &nonterminal = occ.first.first;
    auto &docs = occ.first.second;
    auto &e_poss = occ.second;

    for (std::size_t i = 0; i < docs.size(); ++i) {
      auto pos = dret::FindLastOcc(slp_, nonterminal, docs[i], f_occs_bvs, l_occs_bvs);
      EXPECT_EQ(pos, e_poss[i]) << nonterminal << ":" << docs[i];
    }
  }
}


INSTANTIATE_TEST_SUITE_P(LastOccs,
                         LastOccs_Test,
                         testing::Values(
                             std::make_tuple(2ul,
                                             Rules{{1, 2}, {3, 3}},
                                             Occs{{{1, {1, 2}}, {1, 0}}, {{2, {2, 1}}, {1, 0}}, {{3, {1, 2}}, {1, 2}},
                                                  {{4, {1, 2}}, {3, 4}}}),
                             std::make_tuple(4ul,
                                             Rules{{2, 1}, {3, 5}, {3, 3}, {2, 5}, {4, 6}, {8, 1}, {6, 7}, {11, 6},
                                                   {9, 12}, {13, 10}},
                                             Occs{{{1, {1, 2, 3, 4}}, {1, 0, 0, 0}}, {{2, {1, 2, 3, 4}}, {0, 1, 0, 0}},
                                                  {{3, {1, 2, 3, 4}}, {0, 0, 1, 0}}, {{4, {1, 2, 3, 4}}, {0, 0, 0, 1}},
                                                  {{5, {1, 2}}, {2, 1}}, {{6, {1, 2, 3}}, {3, 2, 1}},
                                                  {{7, {1, 3}}, {0, 2}}, {{10, {1, 2}}, {4, 2}}})
                         )
);
