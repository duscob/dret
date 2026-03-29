//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/27/26.
//


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "dret/doc_list_index_brute.h"

#include "base_test.h"

//~~~~~~~


template <typename TIndex>
class DocListIndexConstructTypedTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(data_);
  }

  const std::string data_ = "MINIMUM\1MINIMAL\1MINIMIZES\1";
};

using DocListIndexConstructTypes = ::testing::Types<dret::DocListIdxBrute<>>;

TYPED_TEST_SUITE(DocListIndexConstructTypedTests, DocListIndexConstructTypes);

TYPED_TEST(DocListIndexConstructTypedTests, construct) {
  auto key_index = "index";
  {
    TypeParam index;
    dret::construct(index, this->config_);
    sdsl::store_to_cache(index, key_index, this->config_, true);
  }

  TypeParam index;
  sdsl::load_from_cache(index, key_index, this->config_, true);
}
