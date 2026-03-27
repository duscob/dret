//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//


#include "base_test.h"
#include "doc_list_index_brute.h"

template <typename TIndex>
class DLIndexTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    data_ = std::make_tuple(String{"MINIMUM\1MINIMAL\1MINIMIZES\1"});

    const auto& data = std::get<0>(data_);
    Init(data);
  }

  std::tuple<String> data_;
};

//~~~~~~~


template <typename TIndex>
class DLIndexConstructTypedTests : public DLIndexTypedTests<TIndex> {};

using DLIndexes = ::testing::Types<dret::DocListIdxBrute<>>;

TYPED_TEST_SUITE(DLIndexConstructTypedTests, DLIndexes);

TYPED_TEST(DLIndexConstructTypedTests, construct) {
  auto key_index = "index";
  {
    TypeParam index;
    dret::construct(index, this->config_);
    sdsl::store_to_cache(index, key_index, this->config_, true);
  }

  TypeParam index;
  sdsl::load_from_cache(index, key_index, this->config_, true);
}
