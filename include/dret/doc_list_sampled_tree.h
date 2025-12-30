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


template <typename TStorage = GenericStorage, typename TAlphabet = Alphabet<>>
class DLSampledTreeScheme : public DocListIndexExtStorage<TStorage, typename TAlphabet::string_type> {
 public:
  using Base = DocListIndexExtStorage<TStorage, typename TAlphabet::string_type>;
  using typename Base::TDocId;
  using typename Base::TPattern;

  DLSampledTreeScheme() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {}
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
