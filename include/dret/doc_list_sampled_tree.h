//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/26/25.
//

#pragma once

#include <filesystem>

#include "sdsl/io.hpp"

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "construct_base.h"
#include "doc_list_index.h"
#include "index_base.h"

namespace dret {

namespace conf {
const std::string KEY_GCDA_CSLP = "gcda_cslp";                // Combined SLP on Document Array
const std::string KEY_GCDA_CSLP_DOCS = "gcda_cslp_docs";      // Chunks for documents
const std::string KEY_GCDA_CSLP_DOCS_C = "gcda_cslp_docs_c";  //

const std::string KEY_GCDA_LSLP = "gcda_lslp";
const std::string KEY_GCDA_LSLP_BASIC = "gcda_lslp_basic";
}  // namespace conf

//~~~~~~~


template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
class DLSampledTreeScheme : public DocListIndexExtStorage<TStorage, typename TAlphabet::string_type> {
 public:
  using Base = DocListIndexExtStorage<TStorage, typename TAlphabet::string_type>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  DLSampledTreeScheme(const TStorage& t_storage) : Base(t_storage), count_idx_(t_storage) {}

  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = count_idx_.Count(t_pattern);

    auto cover = compute_cover_(sp, ep);

    const auto& range = cover.first;
    const auto& nodes = cover.second;

    std::vector<uint32_t> docs;
    docs.reserve(range.first - sp + ep - range.second);

    // TODO Use generic object to copy
    auto add_doc = [&docs](const auto& tt_d) {
      docs.emplace_back(tt_d);
    };

    if (nodes.empty()) {
      get_docs_(sp, ep, add_doc);
    } else {
      get_docs_(sp, range.first, add_doc);
      get_docs_(range.second, ep, add_doc);
    }

    sort(docs.begin(), docs.end());
    docs.erase(unique(docs.begin(), docs.end()), docs.end());

    if (!nodes.empty()) {
      merge_sets_(nodes.begin(), nodes.end(), get_doc_set_, docs);
    }

    for (const auto& doc : docs) {
      t_report(doc);
    }
  }

  void load(Config t_config) override {
    count_idx_.load(t_config);
    compute_cover_.load(t_config);
    get_docs_.load(t_config);
    get_doc_set_.load(t_config);
    merge_sets_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    count_idx_.load(in);
    compute_cover_.load(in);
    get_docs_.load(in);
    get_doc_set_.load(in);
    merge_sets_.load(in);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return count_idx_.serialize(out, child, "count_idx") + compute_cover_.serialize(out, child, "compute_cover")
           + get_docs_.serialize(out, child, "get_docs") + get_doc_set_.serialize(out, child, "get_doc_set")
           + merge_sets_.serialize(out, child, "merge_sets");
  }

 protected:
  TCountIdx count_idx_;
  TComputeCover compute_cover_;
  const TGetDocs get_docs_;
  const TGetDocSet get_doc_set_;
  const TMergeSets merge_sets_;
};

//~~~~~~~


template <typename TStorage = GenericStorage>
class GCDA : public DLSampledTreeScheme<TStorage> {
 public:
  using Base = DLSampledTreeScheme<TStorage>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  GCDA(uint32_t t_block_size, float t_storing_factor) : block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  GCDA() = default;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    // TODO Add implementation
    return 0;
  }

  const uint32_t& block_size() const {
    return block_size_;
  }

  const float& storing_factor() const {
    return storing_factor_;
  }

 protected:
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

//~~~~~~~


void ConstructCombinedSLPOnDA(Config& t_config, uint32_t t_block_size, float t_storing_factor);
void ConstructLightSLPOnDA(Config& t_config);

//~~~~~~~


template <typename TStorage = GenericStorage>
void constructItems(GCDA<TStorage>& t_index, Config& t_config) {
  if (!sdsl::cache_file_exists(sdsl::conf::KEY_SA, t_config)) {
    auto event = sdsl::memory_monitor::event("SA");
    sdsl::construct_sa<8>(t_config);
  }

  if (!cache_file_exists(dret::conf::KEY_DOC_END, t_config)) {
    auto event = sdsl::memory_monitor::event("DocEnds");
    ConstructDocEnd(t_config);
  }

  if (!cache_file_exists(dret::conf::KEY_DA, t_config)) {
    auto event = sdsl::memory_monitor::event("DA");
    ConstructDocArray(t_config);
  }

  if (std::string file_da = cache_file_name(conf::KEY_DA_RAW, t_config);
      !std::filesystem::exists(file_da + ".R") && REPAIR_EXE) {
    auto event = sdsl::memory_monitor::event("DA-Repair");
    std::string cmd = REPAIR_EXE + (" " + file_da);
    std::system(cmd.c_str());
  }

  // Construct Combined SLP of Document Array (Raw)
  if (!sdsl::cache_file_exists(conf::KEY_GCDA_CSLP, t_config)
      || !sdsl::cache_file_exists(conf::KEY_GCDA_CSLP_DOCS, t_config)) {
    auto event = sdsl::memory_monitor::event("DA-CombinedSLP");
    ConstructCombinedSLPOnDA(t_config, t_index.block_size(), t_index.storing_factor());
  }

  // Construct Light SLP of Document Array (Raw)
  if (!sdsl::cache_file_exists(conf::KEY_GCDA_LSLP, t_config)
      || !sdsl::cache_file_exists(conf::KEY_GCDA_LSLP_BASIC, t_config)) {
    auto event = sdsl::memory_monitor::event("DA-LightSLP");
    ConstructLightSLPOnDA(t_config);
  }
}

//~~~~~~~


inline void ConstructCombinedSLPOnDA(Config& t_config, uint32_t t_block_size, float t_storing_factor) {
  auto datafile = cache_file_name(conf::KEY_DA_RAW, t_config);

  grammar::SLP<> slp;
  {
    grammar::RePairReader<true> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(slp);
    re_pair_reader.Read(datafile, slp_wrapper);
  }

  auto cslp = grammar::CombinedSLP<>(slp);
  grammar::Chunks<> cslp_docs;

  grammar::AddSet add_set(cslp_docs);
  cslp.Compute(t_block_size,
               add_set,
               add_set,
               grammar::MustBeSampled<decltype(cslp_docs)>(grammar::AreChildrenTooBig(cslp_docs, t_storing_factor)));

  sdsl::store_to_cache(cslp, conf::KEY_GCDA_CSLP, t_config);

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };

  sdsl::store_to_cache(cslp_docs, conf::KEY_GCDA_CSLP_DOCS, t_config);
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, conf::KEY_GCDA_CSLP_DOCS_C, t_config);
}

void ConstructLightSLPOnDA(Config& t_config) {
  grammar::LightSLP<> lslp;

  {
    // Construct Light SLP on DA
    auto datafile = cache_file_name(conf::KEY_DA_RAW, t_config);

    grammar::SLP<> slp;
    std::vector<std::size_t> compact_seq;
    {
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(slp);

      auto report_compact_seq = [&compact_seq](const auto& _var) {
        compact_seq.emplace_back(_var);
      };

      re_pair_reader.Read(datafile, slp_wrapper, report_compact_seq);
    }

    grammar::CombinedSLP<> cslp;
    sdsl::load_from_cache(cslp, conf::KEY_GCDA_CSLP, t_config);

    lslp.Compute(slp, compact_seq, cslp);
  }

  sdsl::store_to_cache(lslp, conf::KEY_GCDA_LSLP, t_config);

  // Construct Light SLP Basic on DA
  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };
  grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                    grammar::SampledSLP<>,
                    grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>
      lslp_basic(lslp, bit_compress, bit_compress, bit_compress, bit_compress);

  sdsl::store_to_cache(lslp_basic, conf::KEY_GCDA_LSLP_BASIC, t_config);
}

//~~~~~~~

}  // namespace dret
