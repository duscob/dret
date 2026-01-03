//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/9/25.
//

#pragma once

#include "sr-index/sr_idx_generic.h"
#include "sr-index/sr_index.h"

#include "config.h"
#include "construct_base.h"
#include "doc_list_index.h"
#include "index_base.h"

namespace dret {

template <typename TStorage = GenericStorage, typename TBvDocEnds = sdsl::sd_vector<>, uint8_t t_width = 8>
class GetDocBv;

//~~~~~~~


template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TLocateIdx = sri::SrIdxGeneric<sri::SrIndexValidArea<TStorage, TAlphabet>, 8>,
          typename TGetDoc = GetDocBv<TStorage>>
class DocListIdxBrute : public DocListIndexExtStorage<TStorage, typename TAlphabet::string_type> {
 public:
  using Base = DocListIndexExtStorage<TStorage, typename TAlphabet::string_type>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  explicit DocListIdxBrute(const TStorage& t_storage)
      : IndexBaseWithExternalStorage<TStorage>(t_storage), locate_idx_(t_storage), get_doc_(t_storage) {}

  DocListIdxBrute() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto occurrences = locate_idx_.Locate(t_pattern);

    for (const auto& item : occurrences) {
      t_report(get_doc_(item));
    }
  }

  void load(Config t_config) override {
    locate_idx_.load(t_config);
    get_doc_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    locate_idx_.load(in);
    get_doc_.load(in, t_keys);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return locate_idx_.serialize(out, child, "locate_idx") + get_doc_.serialize(out, child, "get_doc");
  }


 protected:
  TLocateIdx locate_idx_;
  TGetDoc get_doc_;
};

//~~~~~~~


template <typename TStorage, typename TLocateIdx, typename TGetDoc>
void constructItems(DocListIdxBrute<TStorage, TLocateIdx, TGetDoc>& t_index, Config& t_config) {
  construct(TLocateIdx{}, t_config);
  construct(TGetDoc{}, t_config);
}

//~~~~~~~


template <typename TStorage, typename TBvDocEnds, uint8_t t_width>
class GetDocBv : public IndexBaseWithExternalStorage<TStorage> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage>;
  using typename Base::size_type;

  GetDocBv(const TStorage& t_storage) : IndexBaseWithExternalStorage<TStorage>(t_storage) {}

  GetDocBv() = default;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TBvDocEnds>(key_, out, child, "doc_end");
    written_bytes += this->template serializeRank<TBvDocEnds>(key_, out, child, "doc_end_rank");

    return written_bytes;
  }

  uint32_t operator()(uint32_t t_position) const {
    return doc_ends_rank_.get(t_position);
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    key_ = t_keys[kDocEnds].get<std::string>();

    this->template loadItem<TBvDocEnds>(key_, t_source, true);
    doc_ends_rank_ = this->template loadBVRank<TBvDocEnds>(key_, t_source, true);
  }

  //~~~~~~~

  std::string key_;
  // std::reference_wrapper<TBvDocEnds> doc_ends_bv_;
  std::reference_wrapper<typename TBvDocEnds::rank_1_type> doc_ends_rank_;
};

//~~~~~~~


template <typename TStorage, typename TBvDocEnds, uint8_t t_width>
void constructItems(GetDocBv<TStorage, TBvDocEnds, t_width>& t_index, Config& t_config) {
  if (!cache_file_exists(sdsl::key_text_trait<t_width>::KEY_TEXT, t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<t_width>(t_config.data_path, t_config);
  }

  if (const auto key = t_config.keys[conf::kDocEnds].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event("DocEnds");
    ConstructDocEnd<t_width, TBvDocEnds>(t_config);
  }
}


}  // namespace dret