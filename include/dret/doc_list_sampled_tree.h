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
}  // namespace conf

//~~~~~~~


template <typename TStorage = GenericStorage>
class DLSampledTreeScheme : public DocListIndex {
 public:
  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {}
};

//~~~~~~~


template <typename TStorage = GenericStorage>
class GCDA : public DLSampledTreeScheme<TStorage> {
 public:
  using Base = DLSampledTreeScheme<TStorage>;
  using typename Base::TDocId;
  using typename Base::TPattern;

  GCDA() = default;

  const uint32_t& block_size() const { return block_size_; }

  const float& storing_factor() const { return storing_factor_; }

 protected:
  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

//~~~~~~~


void ConstructCombinedSLPOnDA(Config& t_config, uint32_t t_block_size, float t_storing_factor);

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
  grammar::CombinedSLP<> cslp;
  if (!sdsl::load_from_cache(cslp, conf::KEY_GCDA_CSLP, t_config)
      || !sdsl::cache_file_exists(conf::KEY_GCDA_CSLP_DOCS, t_config)) {
    auto event = sdsl::memory_monitor::event("DA-CombinedSLP");
    ConstructCombinedSLPOnDA(t_config, t_index.block_size(), t_index.storing_factor());
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
  cslp.Compute(t_block_size, add_set, add_set,
               grammar::MustBeSampled<decltype(cslp_docs)>(grammar::AreChildrenTooBig(cslp_docs, t_storing_factor)));

  sdsl::store_to_cache(cslp, conf::KEY_GCDA_CSLP, t_config);

  auto bit_compress = [](sdsl::int_vector<>& _v) { sdsl::util::bit_compress(_v); };

  sdsl::store_to_cache(cslp_docs, conf::KEY_GCDA_CSLP_DOCS, t_config);
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, conf::KEY_GCDA_CSLP_DOCS_C, t_config);
}

//~~~~~~~

}  // namespace dret
