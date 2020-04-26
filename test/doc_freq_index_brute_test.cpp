//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/25/20.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <doc_freq_index_brute.h>

using Occurrences = std::vector<std::size_t>;
using DocFreqs = std::unordered_map<std::size_t, std::size_t>;

class MockCSA {
 public:
  MockCSA(const Occurrences &_occurrences) : occurrences_{_occurrences} {}

  auto Locate(const std::string &_pattern) const {
    return occurrences_;
  }

 private:
  const Occurrences &occurrences_;
};

class DocFreqIndexBrute_Test : public testing::TestWithParam<std::tuple<Occurrences, DocFreqs>>{};

TEST_P(DocFreqIndexBrute_Test, ListWithFreq) {
  MockCSA csa{std::get<0>(GetParam())};
  auto get_doc = [](const auto &_occ) -> auto {
    return _occ;
  };

  auto idx = dret::MakeDocFreqIndexBrute(csa, get_doc);

  auto real_doc_freqs = idx.ListWithFreq("empty");

  const auto &expected_doc_freqs = std::get<1>(GetParam());

  EXPECT_THAT(real_doc_freqs, testing::UnorderedElementsAreArray(expected_doc_freqs));
}

INSTANTIATE_TEST_SUITE_P(
    DocFreqIndexBrute,
    DocFreqIndexBrute_Test,
    testing::Values(
        std::make_tuple(Occurrences{}, DocFreqs{}),
        std::make_tuple(Occurrences{3}, DocFreqs{{3, 1}}),
        std::make_tuple(Occurrences{3, 3, 1, 3, 1, 6}, DocFreqs{{1, 2}, {3, 3}, {6, 1}}),
        std::make_tuple(Occurrences{3, 3, 1, 3, 1, 6, 321}, DocFreqs{{1, 2}, {3, 3}, {6, 1}, {321, 1}})
    )
);