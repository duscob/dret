//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/3/21.
//

#ifndef DRET_BENCHMARK_DOCUMENT_LISTING_FACTORY_H_
#define DRET_BENCHMARK_DOCUMENT_LISTING_FACTORY_H_

#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/io.hpp>

#include <sr-index/sr_index.h>

#include "../tool/definitions.h"

#include "dret/doc_list_index.h"

using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

template <uint8_t t_width = 8>
class Factory {
 public:
  enum class IndexEnum {
    BRUTE_R_INDEX,
    BRUTE_SR_INDEX,
  };

  struct Config {
    IndexEnum index_t;
  };

  explicit Factory(sri::Config t_config, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : config_{std::move(t_config)} {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_));
    seq_size_ = buf.size();

    r_index_ = std::make_shared<sri::RIndex<ExternalGenericStorage>>(std::ref(storage_));
    r_index_->load(config_);

    sr_index_ = std::make_shared<sri::SrIndexValidArea<ExternalGenericStorage>>(std::ref(storage_), 8);
    sr_index_->load(config_);

    load(doc_endings_);
    load(doc_endings_rank_, [this]() { return TDocEndingsRank(&this->doc_endings_.item); });
    n_doc_ = doc_endings_rank_.item(doc_endings_.item.size());
  }

  std::pair<dret::DocListIndex *, std::size_t> Make(const Config &t_config) {
    return MakeInner(t_config);
  }

  [[nodiscard]] auto SequenceSize() const {
    return seq_size_;
  }

  [[nodiscard]] auto NDocs() const {
    return n_doc_;
  }

 private:

  template<typename T>
  struct Item {
    std::string key;
    bool initialized = false;
    T item;
    std::size_t size_in_bytes = 0;
  };

  template<typename T>
  void load(Item<T> &t_item, const std::string &t_key) {
    if (t_item.initialized) return;

    if (!sdsl::cache_file_exists(t_key, config_))
      std::cerr << "ERROR: File '" << sdsl::cache_file_name(t_key, config_) << "' not exist!!!";

    sdsl::load_from_cache(t_item.item, t_key, config_);
    t_item.initialized = true;
    t_item.size_in_bytes = sdsl::size_in_bytes(t_item.item);
  }

  template<typename T>
  void load(Item<T> &t_item) {
    load(t_item, t_item.key);
  }

  template<typename T, typename TInit>
  void load(Item<T> &t_item, const TInit &t_init) {
    if (t_item.initialized) return;

    t_item.item = t_init();
    t_item.initialized = true;
    t_item.size_in_bytes = sdsl::size_in_bytes(t_item.item);
  }

  std::pair<dret::DocListIndex*, std::size_t> MakeInner(const Config& t_config) {
    dret::DocListIndex* index = nullptr;
    std::size_t index_size = 0;

    switch (t_config.index_t) {
      case IndexEnum::BRUTE_R_INDEX: {
        auto locate = [this](const auto& tt_pattern) { return this->r_index_->Locate(tt_pattern); };

        index = new dret::DocListIndexBrute(locate, doc_endings_rank_.item);
        index_size = sdsl::size_in_bytes(*r_index_) + doc_endings_rank_.size_in_bytes;
        break;
      }

      case IndexEnum::BRUTE_SR_INDEX: {
        auto locate = [this](const auto& tt_pattern) { return this->sr_index_->Locate(tt_pattern); };

        index = new dret::DocListIndexBrute(locate, doc_endings_rank_.item);
        index_size = sdsl::size_in_bytes(*sr_index_) + doc_endings_rank_.size_in_bytes;
        break;
      }
    }

    return std::make_pair(index, index_size);
  }

  sri::Config config_;

  std::size_t seq_size_;

  sri::GenericStorage storage_;

  std::shared_ptr<sri::RIndex<ExternalGenericStorage>> r_index_;
  std::shared_ptr<sri::SrIndexValidArea<ExternalGenericStorage>> sr_index_;

  // Document endings marks
  using TDocEndings = sdsl::sd_vector<>;
  using TDocEndingsRank = TDocEndings::rank_1_type;
  Item<TDocEndings> doc_endings_ = {KEY_DOC_END};
  Item<TDocEndingsRank> doc_endings_rank_;

  // Documents
  std::size_t n_doc_;

};

#endif  //DRET_BENCHMARK_DOCUMENT_LISTING_FACTORY_H_
