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

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
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
          typename TComputeCover = ComputeCover<TStorage>,
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

  DocListIdxGCDA(const TStorage& t_storage, const TCountIdx& t_count_idx, const TComputeCover& t_compute_cover)
      : Base(t_storage, t_count_idx, t_compute_cover) {}

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


template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TComputeCover,
          typename TGetDocs,
          typename TGetDocSet,
          typename TMergeSets>
void construct(DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TComputeCover, TGetDocs, TGetDocSet, TMergeSets>& t_index,
               Config& t_config) {
  if (!cache_file_exists(t_config.keys[conf::kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<TAlphabet::int_width>(t_config);
  }

  auto count_idx = t_index.count_idx;
  construct(count_idx, t_config.data_path, t_config);

  auto compute_cover = t_index.compute_cover;
  construct(compute_cover, t_config);

  t_index.load(t_config);
}

//~~~~~~~


template <typename TStorage, typename TAlphabet, typename TSLP>
class SLPWrapper : public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using typename Base::size_type;

  SLPWrapper(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  SLPWrapper(const TSLP* t_slp, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : slp_(t_slp), block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  SLPWrapper(uint32_t t_block_size, float t_storing_factor)
      : block_size_(t_block_size), storing_factor_(t_storing_factor) {}

  SLPWrapper() = default;

  const TSLP* slp() const {
    return slp_;
  }

  const uint32_t& block_size() const {
    return block_size_;
  }

  const float& storing_factor() const {
    return storing_factor_;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return slp_->serialize(out, child, "slp");
    // return this->template serializeItem<TSLP>(key_, out, child, "slp");
  }

 protected:
  void loadInner(typename Base::TSource& t_source, const JSON& t_keys) override {
    auto key_prefix = std::format("{}-{}_", block_size_, storing_factor_);
    auto key = key_prefix + t_keys[conf::kGCDA][conf::kSLP].get<std::string>();

    slp_ = this->template loadItemPtr<TSLP>(key, t_source, true);
  }

  const TSLP* slp_ = nullptr;

  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

//~~~~~~~


template <typename TStorage, typename TAlphabet, typename TSLP>
class ComputeCover : public SLPWrapper<TStorage, TAlphabet, TSLP> {
 public:
  using Base = SLPWrapper<TStorage, TAlphabet, TSLP>;

  explicit ComputeCover(const TSLP* t_slp, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_slp, t_block_size, t_storing_factor) {}

  explicit ComputeCover(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : Base(t_storage, t_block_size, t_storing_factor) {}

  ComputeCover(uint32_t t_block_size, float t_storing_factor) : Base(t_block_size, t_storing_factor) {}

  ComputeCover() = default;

  auto Compute(std::size_t _bp, std::size_t _ep) const {
    std::vector<std::size_t> nodes;

    auto report = [&nodes](const auto& _value) {
      nodes.emplace_back(_value);
    };

    auto range = grammar::ComputeCoverFromBottom(*this->slp_, _bp, _ep, report);

    return std::make_pair(std::move(range), std::move(nodes));
  }

  auto operator()(std::size_t _sp, std::size_t _ep) const {
    return Compute(_sp, _ep);
  }
};

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TLeavesContainer>
void construct(grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>& t_cslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename TSLP, typename TSampledSLP, typename TChunks>
void construct(grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename TStorage, typename TAlphabet, typename TSLP>
void construct(SLPWrapper<TStorage, TAlphabet, TSLP>& t_slp, Config& t_config) {
  using namespace conf;

  if (const auto key = t_config.keys[kText].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructText<TAlphabet::int_width>(t_config);
  }

  if (const auto key = t_config.keys[kSA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    sdsl::construct_sa<TAlphabet::int_width>(t_config);
  }

  if (const auto key = t_config.keys[kDocEnds].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocEnd<TAlphabet::int_width, sdsl::sd_vector<>>(t_config);
  }

  if (const auto key = t_config.keys[kDA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray(t_config);
  }

  if (const auto key = std::format("{}-{}_", t_slp.block_size(), t_slp.storing_factor())
                       + t_config.keys[kGCDA][kSLP].get<std::string>();
      !sdsl::cache_file_exists<TSLP>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    auto key_da = t_config.keys[conf::kDA].get<std::string>();
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(key_da, t_config);
    TSLP slp = t_slp.slp() ? *t_slp.slp() : TSLP();
    construct(slp, t_config, filepath_da, t_slp.block_size(), t_slp.storing_factor());
  }

  t_slp.load(t_config);
}

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TLeavesContainer>
void construct(grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>& t_cslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;

  // Grammar compress data file using RePair
  if (!std::filesystem::exists(t_datafile + ".R") && REPAIR_EXE) {
    const auto filename = std::filesystem::path(t_datafile).filename().string();
    auto event = sdsl::memory_monitor::event("RePair-" + filename);
    std::string cmd = REPAIR_EXE + (" " + t_datafile);
    std::system(cmd.c_str());
    t_config.file_map[filename + ".R"] = t_datafile + ".R";
    t_config.file_map[filename + ".C"] = t_datafile + ".C";
  }

  const auto key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);

  auto key_slp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>>(key_slp, t_config));

  grammar::SLP<> slp;
  {
    grammar::RePairReader<true> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(slp);
    re_pair_reader.Read(t_datafile, slp_wrapper);
  }

  t_cslp = grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>(slp);
  grammar::Chunks<> cslp_docs;

  grammar::AddSet add_set(cslp_docs);
  t_cslp.Compute(t_block_size,
                 add_set,
                 add_set,
                 grammar::MustBeSampled<decltype(cslp_docs)>(grammar::AreChildrenTooBig(cslp_docs, t_storing_factor)));

  sdsl::store_to_cache(t_cslp, key_slp, t_config, true);

  auto key_docs = key_prefix + t_config.keys[kGCDA][kDocs].get<std::string>();
  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);
}

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TChunks>
void construct(grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;

  std::string key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);

  auto key_lslp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
  grammar::LightSLP<> lslp;

  if (!sdsl::cache_file_exists<decltype(lslp)>(key_lslp, t_config)) {
    // Construct Light SLP on datafile
    auto event = sdsl::memory_monitor::event(sdsl::cache_file_name<grammar::LightSLP<>>(key_lslp, t_config));

    grammar::CombinedSLP<> cslp;
    if (const auto key = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
        !sdsl::cache_file_exists<decltype(cslp)>(key, t_config)) {
      auto event_cslp = sdsl::memory_monitor::event(sdsl::cache_file_name<decltype(cslp)>(key, t_config));
      construct(cslp, t_config, t_datafile, t_block_size, t_storing_factor);
    } else {
      sdsl::load_from_cache(cslp, key, t_config, true);
    }

    grammar::SLP<> slp;
    std::vector<std::size_t> compact_seq;
    {
      grammar::RePairReader<false> re_pair_reader;
      auto slp_wrapper = grammar::BuildSLPWrapper(slp);

      auto report_compact_seq = [&compact_seq](const auto& _var) {
        compact_seq.emplace_back(_var);
      };

      re_pair_reader.Read(t_datafile, slp_wrapper, report_compact_seq);
    }

    lslp.Compute(slp, compact_seq, cslp);
    sdsl::store_to_cache(lslp, key_lslp, t_config, true);
  } else {
    sdsl::load_from_cache(lslp, key_lslp, t_config, true);
  }

  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::LightSLP<TSLP, TSampledSLP, TChunks>>(key_lslp, t_config));

  // Construct Light SLP Basic on DA
  auto bit_compress = [](auto& _v) {
    sdsl::util::bit_compress(_v);
  };
  t_lslp = grammar::LightSLP<TSLP, TSampledSLP, TChunks>(lslp, bit_compress, bit_compress, bit_compress, bit_compress);
  sdsl::store_to_cache(t_lslp, key_lslp, t_config, true);
}

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

//~~~~~~~


class MergeSetsBinaryTreeFunctor {
 public:
  MergeSetsBinaryTreeFunctor() = default;

  template <typename TII, typename TSets, typename TResult>
  void operator()(TII _first, TII _last, const TSets& _sets, TResult& _result) const {
    grammar::MergeSetsBinaryTree(_first, _last, _sets, _result);
  }
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
