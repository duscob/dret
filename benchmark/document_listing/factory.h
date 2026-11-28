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
#include <optional>
#include <string>
#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "sr-index/index_base.h"

#include "../tool/definitions.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/pdl/storage_policy.h"

#include "axes.h"
#include "factories/query.h"


template <uint8_t t_width = 8>
class Factory {
 public:
  // Axis enums are defined at bench::axes scope (axes.h). Aliases below keep
  // every Factory<>::IndexEnum / Factory<>::GetDocEnum / ... reference.
  using IndexEnum        = bench::axes::IndexEnum;
  using PDLVariant       = bench::axes::PDLVariant;
  using PDLStoragePolicy = bench::axes::PDLStoragePolicy;
  using GetDocEnum       = bench::axes::GetDocEnum;
  using GCDASLPVariant   = bench::axes::GCDASLPVariant;
  using BareSLPVariant   = bench::axes::BareSLPVariant;
  using DGCDASLPVariant  = bench::axes::DGCDASLPVariant;
  using RunValuesVariant = bench::axes::RunValuesVariant;
  using PrevDocVariant   = bench::axes::PrevDocVariant;
  using GCDASetsCodec    = bench::axes::GCDASetsCodec;

  // Typed-index aliases for each family live in benchmark/document_listing/
  // factories/{brute,gcda,dgcda,slp_ns,rmq,pdl}.h. The Factory facade
  // dispatches MakeIndex(Config) to the family's query entry point
  // (factories/query.h), compiled in a source file of its own.

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
    // TRunValues axis for the value-rule ILCP-V / CILCP-V cores. Ignored by
    // every other index family. Default DV matches IlcpLikeValueCore's default.
    RunValuesVariant run_values = RunValuesVariant::DV;
    // TPrevDoc axis for the value-rule SADA-V core. Ignored by every other family.
    // Default IV matches SadaValueCore's default (prev_doc values are random
    // SA positions; int_vector with bit_compress is the natural choice).
    PrevDocVariant prev_doc = PrevDocVariant::IV;
    // SA-Phi sampling rate for get_doc=sa_phi_sr (PDL/RMQ × SR-index
    // analogue of RLCSA). Ignored when get_doc != SAPhiSR; r-index variant
    // (sa_phi_r) has no sampling knob.
    std::size_t sa_sampling = 0;
    // Block size of PDL's GCDA backend (get_doc=slp, light), the grammar that
    // answers the ranges PDL's sampled nodes do not cover. 1024 is the value it
    // was pinned at; its storing factor stays 32 (it only shapes the document
    // lists, which the backend does not load). Ignored by every other family.
    uint32_t pdl_backend_block_size = 1024;
    // How GCDA / DGCDA store their sampled nodes' document sets (Re-Pair or
    // plain). Ignored by every other family.
    GCDASetsCodec gcda_sets = GCDASetsCodec::RP;
    // SLP_NS over a sampled GCDA grammar (gcda_slp at block_size,
    // storing_factor) instead of a bare one (bare_slp). Ignored by every other
    // family.
    bool slp_ns_sampled = false;

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
      if (pdl_backend_block_size != t_c.pdl_backend_block_size)
        return pdl_backend_block_size < t_c.pdl_backend_block_size;
      if (gcda_sets != t_c.gcda_sets)
        return gcda_sets < t_c.gcda_sets;
      if (slp_ns_sampled != t_c.slp_ns_sampled)
        return slp_ns_sampled < t_c.slp_ns_sampled;
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

    // The r-index CountOccurrences reads is NOT loaded eagerly: it backs only
    // the per-pattern occurrence counts and requires the full r-index sample /
    // mark cache files (typed bwt_run_first_text_pos_<HASH> etc.). Loading it
    // in the ctor would force every Factory user -- including PDL-only or
    // GCDA-only benchmark runs that never touch the r-index -- to pre-build
    // those artifacts. Deferred to the first CountOccurrences call.

    load(doc_endings_);
    load(doc_endings_rank_, [this]() {
      return TDocEndingsRank(&this->doc_endings_.item);
    });
    n_doc_ = doc_endings_rank_.item(doc_endings_.item.size());
  }

  // Occurrences of t_pattern in the text: the size of its suffix-array range in the
  // brute baseline's r-index (half-open -- RIndexCount treats start >= end as empty).
  // Loads the r-index on first use, so it needs the same cache files as BruteRI.
  std::size_t CountOccurrences(const std::string& t_pattern) {
    if (!count_occurrences_)
      count_occurrences_ = bench::factories::query::OccurrenceCounter(std::ref(storage_), config_);
    return count_occurrences_(t_pattern);
  }

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
  // parsing) can use the same mapping.
  static dret::pdl::StoragePolicy toPDLStoragePolicy(PDLStoragePolicy t_p) {
    switch (t_p) {
      case PDLStoragePolicy::StoreAllInternal:
        return dret::pdl::StoragePolicy::StoreAllInternal;
      case PDLStoragePolicy::LeavesOnly:
        return dret::pdl::StoragePolicy::LeavesOnly;
      case PDLStoragePolicy::ListWeighted:
        return dret::pdl::StoragePolicy::ListWeighted;
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

    namespace q = bench::factories::query;
    const q::RMQArgs rmq{.block_size = t_config.block_size,
                         .storing_factor = t_config.storing_factor,
                         .get_doc = t_config.get_doc,
                         .gcda_slp = t_config.gcda_slp,
                         .bare_slp = t_config.bare_slp,
                         .dgcda_slp = t_config.dgcda_slp,
                         .run_values = t_config.run_values,
                         .prev_doc = t_config.prev_doc,
                         .sa_sampling = t_config.sa_sampling};
    const q::PDLArgs pdl{.block_size = t_config.block_size,
                         .storing_factor = t_config.storing_factor,
                         .get_doc = t_config.get_doc,
                         .policy = t_config.pdl_storage_policy,
                         .gcda_slp = t_config.gcda_slp,
                         .bare_slp = t_config.bare_slp,
                         .dgcda_slp = t_config.dgcda_slp,
                         .sa_sampling = t_config.sa_sampling,
                         .backend_block_size = t_config.pdl_backend_block_size};
    const auto storage = std::ref(storage_);

    Index index;
    q::Built built;
    switch (t_config.index_t) {
      // Brute baselines bypass auto-build: the underlying sri::RIndex /
      // sri::SrIndexValidArea caches come from a construct run, and
      // DocListIdxBrute's own construct() would rebuild them from scratch.
      case IndexEnum::BRUTE_R_INDEX:
        built = q::BruteR(storage, config_);
        break;
      case IndexEnum::BRUTE_SR_INDEX:
        built = q::BruteSr(storage, config_, t_config.sampling_size);
        break;
      case IndexEnum::BRUTE_DA:
        built = q::BruteDA(storage, config_);
        break;
      case IndexEnum::COUNT_RI:
        built = q::CountRI(storage, config_);
        break;
      case IndexEnum::BRUTE_RLCSA:
        built = q::BruteRLCSA(storage, config_);
        break;

      case IndexEnum::GCDA:
        built = q::GCDA(
            storage, config_, t_config.block_size, t_config.storing_factor, t_config.gcda_slp, t_config.gcda_sets);
        break;
      case IndexEnum::GCDA_BACKEND:
        built = q::GCDABackend(
            storage, config_, t_config.block_size, t_config.storing_factor, t_config.get_doc, t_config.gcda_sets);
        break;
      case IndexEnum::DGCDA:
        built = q::DGCDA(
            storage, config_, t_config.block_size, t_config.storing_factor, t_config.dgcda_slp, t_config.gcda_sets);
        break;
      case IndexEnum::SLP_NS:
        built = q::SLPNS(storage,
                         config_,
                         t_config.bare_slp,
                         t_config.block_size,
                         t_config.slp_ns_sampled ? std::optional(t_config.gcda_slp) : std::nullopt,
                         t_config.storing_factor);
        break;

      // Each core reads the axes it has: the published ones neither run_values nor prev_doc,
      // SADA-V prev_doc, ILCP-V and CILCP-V run_values.
      case IndexEnum::SADA:
        built = q::Sada(storage, config_, rmq);
        break;
      case IndexEnum::ILCP:
        built = q::Ilcp(storage, config_, rmq);
        break;
      case IndexEnum::CILCP:
        built = q::Cilcp(storage, config_, rmq);
        break;
      case IndexEnum::SADA_V:
        built = q::SadaValue(storage, config_, rmq);
        break;
      case IndexEnum::ILCP_V:
        built = q::IlcpValue(storage, config_, rmq);
        break;
      case IndexEnum::CILCP_V:
        built = q::CilcpValue(storage, config_, rmq);
        break;

      case IndexEnum::PDL:
        switch (t_config.pdl_variant) {
          case PDLVariant::RP:
            built = q::PDLRP(storage, config_, pdl);
            break;
          case PDLVariant::BC:
            built = q::PDLBC(storage, config_, pdl);
            break;
          case PDLVariant::Plain:
          default:
            built = q::PDLPlain(storage, config_, pdl);
            break;
        }
        break;
    }
    index = {built.first, built.second};

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

  dret::Config config_;

  std::size_t seq_size_;

  sri::GenericStorage storage_;

  std::map<Config, Index> indexes_;

  // Pattern occurrences, through an r-index loaded on first use (factories/brute.cpp).
  std::function<std::size_t(const std::string&)> count_occurrences_;

  // Document endings marks
  using TDocEndings = sdsl::sd_vector<>;
  using TDocEndingsRank = TDocEndings::rank_1_type;
  Item<TDocEndings> doc_endings_ = {KEY_DOC_END};
  Item<TDocEndingsRank> doc_endings_rank_;

  // Documents
  std::size_t n_doc_;
};
