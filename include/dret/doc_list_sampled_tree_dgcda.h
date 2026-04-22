//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/21/26.
//

#pragma once

#include <sdsl/construct_sa.hpp>

#include <grammar/re_pair.h>

#include "construct_base.h"
#include "differential_light_slp.h"
#include "doc_list_sampled_tree_gcda.h"

namespace dret {

namespace dgcda {

template <typename TStorage   = GenericStorage,
          typename TAlphabet  = Alphabet<>,
          typename TCountIdx  = sri::SrIdxGeneric<sri::SrIndexValidArea<TStorage, TAlphabet>, 16>,
          typename TSLP       = DifferentialLightSLP<>,
          typename TSLPSets   = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                                   true,
                                                   grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
class DocListIdxDGCDA
    : public gcda::DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TSLP, TSLPSets, TMergeSets> {
  using Base = gcda::DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TSLP, TSLPSets, TMergeSets>;
 public:
  using Base::Base;

 protected:
  void loadInner(typename Base::StorageBase::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit(
        [this](auto&& tt_source) {
          this->count_idx_.load(tt_source.get());
        },
        t_source);

    auto key_prefix = std::format("{}-{}_", this->block_size_, this->storing_factor_);
    this->slp_ = this->template loadItemPtr<TSLP>(
        key_prefix + t_keys[kDGCDA][kSLP].get<std::string>(), t_source, true);
    this->slp_sets_ = this->template loadItemPtr<TSLPSets>(
        key_prefix + t_keys[kDGCDA][kDocs].get<std::string>(), t_source, true);
  }
};

//~~~~~~~


template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TSLP,
          typename TSLPSets,
          typename TMergeSets>
void construct(DocListIdxDGCDA<TStorage, TAlphabet, TCountIdx, TSLP, TSLPSets, TMergeSets>& t_index,
               Config& t_config) {
  using namespace conf;

  if (!cache_file_exists(t_config.keys[kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
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

  const auto key_prefix = std::format("{}-{}_", t_index.block_size(), t_index.storing_factor());

  // Build DifferentialLightSLP (also stores plain Chunks<> under kDGCDA/kDocs)
  if (const auto key_dslp = key_prefix + t_config.keys[kDGCDA][kSLP].get<std::string>();
      !sdsl::cache_file_exists<TSLP>(key_dslp, t_config)) {
    auto event = sdsl::memory_monitor::event(key_dslp);
    TSLP dslp;
    construct(dslp, t_config, t_index.block_size(), t_index.storing_factor());
  }

  // Build GCChunks from the stored plain Chunks<> under kDGCDA/kDocs
  const auto key_docs = key_prefix + t_config.keys[kDGCDA][kDocs].get<std::string>();
  if (!sdsl::cache_file_exists<TSLPSets>(key_docs, t_config)) {
    auto event = sdsl::memory_monitor::event(key_docs);

    grammar::Chunks<> cslp_docs;
    sdsl::load_from_cache(cslp_docs, key_docs, t_config, true);

    auto bit_compress = [](sdsl::int_vector<>& v) { sdsl::util::bit_compress(v); };
    const auto& objs = cslp_docs.GetObjects();

    grammar::GCChunks<grammar::SLP<>> gc_slp;
    grammar::RePairEncoder<false> encoder_nslp;
    gc_slp.Compute(objs.begin(), objs.end(), cslp_docs, encoder_nslp);
    sdsl::store_to_cache(gc_slp, key_docs, t_config, true);

    TSLPSets slp_sets(gc_slp, bit_compress, bit_compress, bit_compress, bit_compress);
    sdsl::store_to_cache(slp_sets, key_docs, t_config, true);
  }

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  t_index.load(t_config);
}

}  // namespace dgcda

}  // namespace dret
