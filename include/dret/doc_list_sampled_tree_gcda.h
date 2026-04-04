//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/21/26.
//

#pragma once

#include <sdsl/construct_sa.hpp>

#include <grammar/slp_helper.h>

#include "doc_list_sampled_tree.h"
#include "slp_tools.h"

namespace dret {

namespace conf {
const std::string KEY_GCDA_CSLP = "gcda_cslp";                // Combined SLP on Document Array
const std::string KEY_GCDA_CSLP_DOCS = "gcda_cslp_docs";      // Chunks for documents
const std::string KEY_GCDA_CSLP_DOCS_C = "gcda_cslp_docs_c";  //

const std::string KEY_GCDA_LSLP = "gcda_lslp";
const std::string KEY_GCDA_LSLP_BASIC = "gcda_lslp_basic";
}  // namespace conf

//~~~~~~~


namespace gcda {

template <typename TAlphabet = Alphabet<>,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>>
class ComputeCover;

template <typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>>
class GetDocs;

template <typename TSLP = grammar::BasicSLP<sdsl::int_vector<>>,
          bool kExpand = true,
          typename TChunks = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>
class GetDocSet;

class MergeSetsBinaryTreeFunctor;

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::SrIdxGeneric<sri::SrIndexValidArea<TStorage, TAlphabet>, 16>,
          typename TComputeCover = ComputeCover<>,
          typename TGetDocs = GetDocs<>,
          typename TGetDocSet = GetDocSet<>,
          typename TMergeSets = MergeSetsBinaryTreeFunctor>
class DocListIdxGCDA
    : public DLSampledTreeScheme<TStorage, TAlphabet, TCountIdx, TComputeCover, TGetDocs, TGetDocSet, TMergeSets> {
 public:
  using Base = DLSampledTreeScheme<TStorage, TAlphabet, TCountIdx, TComputeCover, TGetDocs, TGetDocSet, TMergeSets>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  explicit DocListIdxGCDA(const TStorage& t_storage) : Base(t_storage) {}

  DocListIdxGCDA(uint32_t t_block_size, float t_storing_factor)
      : block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  DocListIdxGCDA() = default;

  // size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
  //   // TODO Add implementation
  //   return 0;
  // }

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


template <typename TAlphabet, typename TSLP>
class ComputeCover {
 public:
  explicit ComputeCover(const TSLP& _slp) : slp_(_slp) {}

  ComputeCover() = default;

  auto Compute(std::size_t _bp, std::size_t _ep) const {
    std::vector<std::size_t> nodes;

    auto report = [&nodes](const auto& _value) {
      nodes.emplace_back(_value);
    };

    auto range = grammar::ComputeCoverFromBottom(slp_, _bp, _ep, report);

    return std::make_pair(std::move(range), std::move(nodes));
  }

  auto operator()(std::size_t _sp, std::size_t _ep) const {
    return Compute(_sp, _ep);
  }

 protected:
  TSLP slp_;
};


//~~~~~~~


template <typename TSLP>
class GetDocs {
 public:
  GetDocs(const TSLP& _slp) : slp_(_slp) {}

  GetDocs() = default;

  template <typename Report>
  void operator()(std::size_t _bp, std::size_t _ep, Report& _report) const {
    ExpandSLP(slp_, _bp, _ep, _report);
  }

 private:
  TSLP slp_;
};

//~~~~~~~


template <typename TSLP, bool kExpand, typename TChunks>
class GetDocSet : public grammar::GCChunks<TSLP, kExpand, TChunks> {
 public:
  using Base = grammar::GCChunks<TSLP, kExpand, TChunks>;

  GetDocSet() = default;
};

}  // namespace gcda

//~~~~~~~


void ConstructCombinedSLPOnDA(Config& t_config, uint32_t t_block_size, float t_storing_factor);
void ConstructLightSLPOnDA(Config& t_config);

//~~~~~~~


template <typename TStorage = GenericStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
void constructItems(
    gcda::DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TComputeCover, TGetDocs, TGetDocSet, TMergeSets>& t_index,
    Config& t_config) {
  if (!cache_file_exists(t_config.keys[conf::kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<TAlphabet::int_width>(t_config);
  }

  if (const auto key = t_config.keys[conf::kSA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    sdsl::construct_sa<8>(t_config);
  }

  if (const auto key = t_config.keys[conf::kDocEnds].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocEnd(t_config);
  }

  if (const auto key = t_config.keys[conf::kDA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
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
