//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/3/21.
//

#pragma once


#include <any>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "../tool/definitions.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"

#include "axes.h"
#include "factories/brute.h"
#include "factories/dgcda.h"
#include "factories/gcda.h"
#include "factories/pdl.h"
#include "factories/rmq.h"
#include "factories/slp_ns.h"


using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

template <uint8_t t_width = 8>
class Factory {
 public:
  // Axis enums are defined at bench::axes scope (axes.h) so bm_query_doc_list
  // and bm_build_items can share the same types. Aliases below preserve every
  // existing Factory<>::IndexEnum / Factory<>::GetDocEnum / ... reference.
  using IndexEnum        = bench::axes::IndexEnum;
  using PDLVariant       = bench::axes::PDLVariant;
  using PDLStoragePolicy = bench::axes::PDLStoragePolicy;
  using GetDocEnum       = bench::axes::GetDocEnum;
  using GCDASLPVariant   = bench::axes::GCDASLPVariant;
  using BareSLPVariant   = bench::axes::BareSLPVariant;
  using DGCDASLPVariant  = bench::axes::DGCDASLPVariant;
  using RunValuesVariant = bench::axes::RunValuesVariant;
  using PrevDocVariant   = bench::axes::PrevDocVariant;

  // Typed-index aliases for each family now live in benchmark/document_listing/
  // factories/{brute,gcda,dgcda,slp_ns,rmq,pdl}.h. The Factory facade simply
  // dispatches MakeIndex(Config) to the right family's Make() function.

  struct Config {
    IndexEnum index_t;
    std::size_t sampling_size = 0;
    uint32_t block_size = 512;
    float storing_factor = 4;
    GetDocEnum get_doc = GetDocEnum::DA;
    GCDASLPVariant gcda_slp = GCDASLPVariant::Light;
    BareSLPVariant bare_slp = BareSLPVariant::IV;
    PDLVariant pdl_variant = PDLVariant::Plain;
    PDLStoragePolicy pdl_storage_policy = PDLStoragePolicy::OccurrenceWeighted;
    // Appended at the end of the struct so older 9-arg positional Config{}
    // initialisers in the benchmark binaries keep landing in the correct
    // fields. New code should prefer designated init.
    DGCDASLPVariant dgcda_slp = DGCDASLPVariant::Default;
    // TRunValues axis for the published ILCP / CILCP cores. Ignored by
    // every other index family. Default DV matches IlcpLikeFullCore's default.
    RunValuesVariant run_values = RunValuesVariant::DV;
    // TPrevDoc axis for the published SADA core. Ignored by every other family.
    // Default IV matches SadaCore's default (prev_doc values are random
    // SA positions; int_vector with bit_compress is the natural choice).
    PrevDocVariant prev_doc = PrevDocVariant::IV;
    // SA-Phi sampling rate for get_doc=sa_phi_sr (PDL/RMQ × SR-index
    // analogue of RLCSA). Ignored when get_doc != SAPhiSR; r-index variant
    // (sa_phi_r) has no sampling knob.
    std::size_t sa_sampling = 0;

    bool operator<(const Config& t_c) const {
      if (index_t != t_c.index_t)
        return index_t < t_c.index_t;
      if (sampling_size != t_c.sampling_size)
        return sampling_size < t_c.sampling_size;
      if (block_size != t_c.block_size)
        return block_size < t_c.block_size;
      if (storing_factor != t_c.storing_factor)
        return storing_factor < t_c.storing_factor;
      if (get_doc != t_c.get_doc)
        return get_doc < t_c.get_doc;
      if (gcda_slp != t_c.gcda_slp)
        return gcda_slp < t_c.gcda_slp;
      if (bare_slp != t_c.bare_slp)
        return bare_slp < t_c.bare_slp;
      if (pdl_variant != t_c.pdl_variant)
        return pdl_variant < t_c.pdl_variant;
      if (pdl_storage_policy != t_c.pdl_storage_policy)
        return pdl_storage_policy < t_c.pdl_storage_policy;
      if (dgcda_slp != t_c.dgcda_slp)
        return dgcda_slp < t_c.dgcda_slp;
      if (run_values != t_c.run_values)
        return run_values < t_c.run_values;
      if (prev_doc != t_c.prev_doc)
        return prev_doc < t_c.prev_doc;
      return sa_sampling < t_c.sa_sampling;
    }
  };

  explicit Factory(dret::Config t_config, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : config_{std::move(t_config)} {
    // Sequence length, read straight from the BWT's int_vector header.
    //
    // This used to be `int_vector_buffer<t_width> buf(...); seq_size_ = buf.size();`
    // which REWROTE the file on every query run. sdsl::int_vector_buffer's
    // constructor forces the output stream open regardless of the requested
    // mode --- `m_ofile.open(m_filename, mode | std::ios::out | std::ios::binary)`
    // --- and its destructor calls close(), which unconditionally does
    // write_block() and rewrites the header. So merely asking a read-only-looking
    // buffer for its size touched bwt_data.sdsl every time a Factory was built.
    //
    // Harmless in content (the BWT is deterministic, and the 8 collections whose
    // copies happened to be read-only completed identically with the write
    // failing), but it dirtied mtimes, which is what makes "did the query stage
    // build anything?" checks unreliable --- they compare filenames, not mtimes.
    //
    // Mirrors int_vector_buffer's own arithmetic: read_header yields the length
    // in BITS and the stored width, and m_size = bits / width.
    {
      const auto bwt_file = sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_);
      std::ifstream bwt_in(bwt_file, std::ios::binary);
      if (!bwt_in) {
        throw std::runtime_error("Factory: cannot open BWT file '" + bwt_file + "'");
      }
      uint64_t bits = 0;
      uint8_t stored_width = 0;
      sdsl::int_vector<0>::read_header(bits, stored_width, bwt_in);
      if (stored_width == 0) {
        throw std::runtime_error("Factory: BWT '" + bwt_file + "' reports width 0");
      }
      seq_size_ = bits / stored_width;
    }

    // r_index_ / sr_index_ are NOT loaded eagerly: they back only the
    // BRUTE_R_INDEX / BRUTE_SR_INDEX paths in MakeInner and require
    // the full sr-index sample / mark cache files (typed
    // bwt_run_first_text_pos_<HASH> etc.). Loading them in the ctor
    // would force every Factory user — including PDL-only or
    // GCDA-only benchmark runs that never touch the r-index — to
    // pre-build those artifacts. Deferred to rIndex() / srIndex().

    load(doc_endings_);
    load(doc_endings_rank_, [this]() {
      return TDocEndingsRank(&this->doc_endings_.item);
    });
    n_doc_ = doc_endings_rank_.item(doc_endings_.item.size());
  }

 private:
  // Lazy accessors for the brute baselines' shared r-index / sr-index
  // instances. First call constructs + loads from the cache; later
  // calls return the cached shared_ptr.
  sri::RIndex<ExternalGenericStorage>& rIndex() {
    if (!r_index_) {
      r_index_ = std::make_shared<sri::RIndex<ExternalGenericStorage>>(std::ref(storage_));
      r_index_->load(config_);
    }
    return *r_index_;
  }

  sri::SrIndexValidArea<ExternalGenericStorage>& srIndex() {
    if (!sr_index_) {
      sr_index_ = std::make_shared<sri::SrIndexValidArea<ExternalGenericStorage>>(std::ref(storage_), 8);
      sr_index_->load(config_);
    }
    return *sr_index_;
  }

 public:

  std::pair<dret::DocListIndex<>*, std::size_t> Make(const Config& t_config) {
    auto index = MakeIndex(t_config);
    return {index.idx.get(), index.size};
  }

  [[nodiscard]] auto SequenceSize() const {
    return seq_size_;
  }

  [[nodiscard]] auto NDocs() const {
    return n_doc_;
  }

  struct Index {
    std::shared_ptr<dret::DocListIndex<>> idx;
    std::size_t size = 0;
  };

  // Map the factory-side PDLStoragePolicy enum to the library-side
  // dret::pdl::StoragePolicy. Kept as a static so callers (e.g. CLI
  // parsing in Task 30) can use the same mapping.
  static dret::pdl::StoragePolicy toPDLStoragePolicy(PDLStoragePolicy t_p) {
    switch (t_p) {
      case PDLStoragePolicy::StoreAllInternal:
        return dret::pdl::StoragePolicy::StoreAllInternal;
      case PDLStoragePolicy::LeavesOnly:
        return dret::pdl::StoragePolicy::LeavesOnly;
      case PDLStoragePolicy::OccurrenceWeighted:
      default:
        return dret::pdl::StoragePolicy::OccurrenceWeighted;
    }
  }

  Index MakeIndex(const Config& t_config) {
    auto it = indexes_.find(t_config);
    if (it != indexes_.end()) {
      return it->second;
    }

    Index index;
    switch (t_config.index_t) {
      case IndexEnum::BRUTE_R_INDEX: {
        // Brute baselines bypass auto-build: the underlying sri::RIndex /
        // sri::SrIndexValidArea caches are managed externally via the
        // factory's lazy rIndex() / srIndex(), and DocListIdxBrute's own
        // construct() would try to rebuild them from scratch.
        auto idx = std::make_shared<bench::factories::brute::IdxR<ExternalGenericStorage>>(
            std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::BRUTE_SR_INDEX: {
        auto idx = std::make_shared<bench::factories::brute::IdxSr<ExternalGenericStorage>>(
            std::ref(storage_),
            sri::SrIndexValidArea<ExternalGenericStorage>(storage_, t_config.sampling_size));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::GCDA: {
        auto [idx, size] = bench::factories::gcda::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor, t_config.gcda_slp);
        index = {idx, size};
        break;
      }

      case IndexEnum::DGCDA: {
        auto [idx, size] = bench::factories::dgcda::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor, t_config.dgcda_slp);
        index = {idx, size};
        break;
      }

      case IndexEnum::SADA: {
        auto [idx, size] = bench::factories::rmq::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            bench::factories::rmq::CoreKind::SADA_L,
            t_config.get_doc, t_config.gcda_slp, t_config.bare_slp, t_config.dgcda_slp,
            bench::axes::RunValuesVariant::DV, bench::axes::PrevDocVariant::IV,
            t_config.sa_sampling);
        index = {idx, size};
        break;
      }

      case IndexEnum::ILCP: {
        auto [idx, size] = bench::factories::rmq::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            bench::factories::rmq::CoreKind::ILCP_L,
            t_config.get_doc, t_config.gcda_slp, t_config.bare_slp, t_config.dgcda_slp,
            bench::axes::RunValuesVariant::DV, bench::axes::PrevDocVariant::IV,
            t_config.sa_sampling);
        index = {idx, size};
        break;
      }

      case IndexEnum::PDL: {
        auto [idx, size] = bench::factories::pdl::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            t_config.pdl_variant, t_config.get_doc, t_config.pdl_storage_policy,
            t_config.gcda_slp, t_config.bare_slp, t_config.dgcda_slp,
            t_config.sa_sampling);
        index = {idx, size};
        break;
      }

      case IndexEnum::SLP_NS: {
        auto [idx, size] = bench::factories::slp_ns::Make(
            std::ref(storage_), config_, t_config.bare_slp);
        index = {idx, size};
        break;
      }

      case IndexEnum::CILCP: {
        auto [idx, size] = bench::factories::rmq::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            bench::factories::rmq::CoreKind::CILCP_L,
            t_config.get_doc, t_config.gcda_slp, t_config.bare_slp, t_config.dgcda_slp,
            bench::axes::RunValuesVariant::DV, bench::axes::PrevDocVariant::IV,
            t_config.sa_sampling);
        index = {idx, size};
        break;
      }

      case IndexEnum::SADA_S: {
        auto [idx, size] = bench::factories::rmq::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            bench::factories::rmq::CoreKind::SADA,
            t_config.get_doc, t_config.gcda_slp, t_config.bare_slp,
            t_config.dgcda_slp, t_config.run_values, t_config.prev_doc,
            t_config.sa_sampling);
        index = {idx, size};
        break;
      }

      case IndexEnum::ILCP_S: {
        auto [idx, size] = bench::factories::rmq::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            bench::factories::rmq::CoreKind::ILCP,
            t_config.get_doc, t_config.gcda_slp, t_config.bare_slp,
            t_config.dgcda_slp, t_config.run_values,
            bench::axes::PrevDocVariant::IV, t_config.sa_sampling);
        index = {idx, size};
        break;
      }

      case IndexEnum::CILCP_S: {
        auto [idx, size] = bench::factories::rmq::Make(
            std::ref(storage_), config_,
            t_config.block_size, t_config.storing_factor,
            bench::factories::rmq::CoreKind::CILCP,
            t_config.get_doc, t_config.gcda_slp, t_config.bare_slp,
            t_config.dgcda_slp, t_config.run_values,
            bench::axes::PrevDocVariant::IV, t_config.sa_sampling);
        index = {idx, size};
        break;
      }
    }

    if (index.idx) {
      indexes_[t_config] = index;
    }

    return index;
  }


 private:
  template <typename T>
  struct Item {
    std::string key;
    bool initialized = false;
    T item;
    std::size_t size_in_bytes = 0;
  };

  template <typename T>
  void load(Item<T>& t_item, const std::string& t_key) {
    if (t_item.initialized)
      return;

    if (!sdsl::cache_file_exists(t_key, config_))
      std::cerr << "ERROR: File '" << sdsl::cache_file_name(t_key, config_) << "' not exist!!!";

    sdsl::load_from_cache(t_item.item, t_key, config_);
    t_item.initialized = true;
    t_item.size_in_bytes = sdsl::size_in_bytes(t_item.item);
  }

  template <typename T>
  void load(Item<T>& t_item) {
    load(t_item, t_item.key);
  }

  template <typename T, typename TInit>
  void load(Item<T>& t_item, const TInit& t_init) {
    if (t_item.initialized)
      return;

    t_item.item = t_init();
    t_item.initialized = true;
    t_item.size_in_bytes = sdsl::size_in_bytes(t_item.item);
  }

  std::pair<dret::DocListIndex<>*, std::size_t> MakeInner(const Config& t_config) {
    dret::DocListIndex<>* index = nullptr;
    std::size_t index_size = 0;

    switch (t_config.index_t) {
      case IndexEnum::BRUTE_R_INDEX: {
        rIndex();  // trigger lazy load of r_index_ before capturing this.
        auto locate = [this](const auto& tt_pattern) {
          return this->r_index_->Locate(tt_pattern);
        };

        index = new dret::DocListIndexBrute(locate, doc_endings_rank_.item);
        index_size = sdsl::size_in_bytes(*r_index_) + doc_endings_rank_.size_in_bytes;
        break;
      }

      case IndexEnum::BRUTE_SR_INDEX: {
        srIndex();  // trigger lazy load of sr_index_ before capturing this.
        auto locate = [this](const auto& tt_pattern) {
          return this->sr_index_->Locate(tt_pattern);
        };

        index = new dret::DocListIndexBrute(locate, doc_endings_rank_.item);
        index_size = sdsl::size_in_bytes(*sr_index_) + doc_endings_rank_.size_in_bytes;
        break;
      }
    }

    return std::make_pair(index, index_size);
  }

  dret::Config config_;

  std::size_t seq_size_;

  sri::GenericStorage storage_;

  std::map<Config, Index> indexes_;

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
