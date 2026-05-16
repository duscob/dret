//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/3/21.
//

#pragma once


#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/vlc_vector.hpp>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "../tool/definitions.h"

#include "dret/slp/basic_slp_span_length.h"
#include "dret/config.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/doc_list/doc_list_slp.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/doc_list/doc_list_brute.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/doc_list/doc_list_gcda.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/doc_list/doc_list_pdl.h"
#include "dret/pdl/get_docs.h"
#include "dret/pdl/storage_policy.h"

#include "axes.h"


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

  // DGCDA variants differ only in the DifferentialLightSLP's inner TSLP type
  // (the span-length strategy). Everything downstream — SampledSLP, Chunks,
  // GCChunks, count index — is unchanged. Cache files are type-hashed by sdsl,
  // so the variants do not collide on disk.
  using DGCDASLP_OTF = dret::DifferentialLightSLP<
      dret::BasicSLPOnTheFlySpanLength<grammar::BasicSLP<>>>;
  using DGCDASLP_CRL = dret::DifferentialLightSLP<
      dret::BasicSLPCachedRootSpanLengths<grammar::BasicSLP<>>>;

  // Vary the three non-monotonic int-vector fields (roots_, span_sums_, samples_).
  // sample_roots_pos_ keeps its class default sdsl::enc_vector<> (strictly monotonic).
  using DGCDASLP_EV = dret::DifferentialLightSLP<grammar::SLP<>,
                                                  grammar::SampledSLP<>,
                                                  sdsl::enc_vector<>,
                                                  sdsl::enc_vector<>,
                                                  sdsl::enc_vector<>>;
  using DGCDASLP_DV = dret::DifferentialLightSLP<grammar::SLP<>,
                                                  grammar::SampledSLP<>,
                                                  sdsl::dac_vector<>,
                                                  sdsl::dac_vector<>,
                                                  sdsl::dac_vector<>>;
  using DGCDASLP_VV = dret::DifferentialLightSLP<grammar::SLP<>,
                                                  grammar::SampledSLP<>,
                                                  sdsl::vlc_vector<>,
                                                  sdsl::vlc_vector<>,
                                                  sdsl::vlc_vector<>>;

  template <typename TSLP>
  using DGCDAVariant = dret::dgcda::DocListIdxDGCDA<
      ExternalGenericStorage,
      dret::Alphabet<>,
      sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
      TSLP>;

  // GCDA TSLP variants (Phase A — Compact-BP / Compact-LOUDS; Phase B — CSLP).
  using GCDASLP_CompactBP    = grammar::CompactBPSLP<>;
  using GCDASLP_CompactLOUDS = grammar::CompactLOUDSSLP<>;
  using GCDASLP_CSLP         = grammar::CombinedSLPWithUnitCover<>;

  template <typename TSLP>
  using GCDAVariant = dret::gcda::DocListIdxGCDA<
      ExternalGenericStorage,
      dret::Alphabet<>,
      sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
      TSLP>;

  using TCountIdx = sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>;
  using GetDocSLP = dret::rmq::GetDocSLP<ExternalGenericStorage>;
  using GetDocSLP_NS = dret::rmq::GetDocSLP_NS<ExternalGenericStorage>;
  using GetDocDSLP = dret::rmq::GetDocDSLP<ExternalGenericStorage>;

  // Parametrized GetDocSLP for compact-grammar TSLP variants. The default
  // GetDocSLP defaults its TSLP to grammar::LightSLP<...> (matching default
  // GCDA), so cache sharing is automatic; here we deliberately use different
  // TSLPs so SADA-SLP / ILCP-SLP / CILCP-SLP can reuse the SLP cache files
  // produced for the corresponding GCDA-CompactBP / GCDA-CompactLOUDS build.
  using GetDocSLP_CompactBP    = dret::rmq::GetDocSLP<ExternalGenericStorage,
                                                       dret::Alphabet<>::int_width,
                                                       GCDASLP_CompactBP>;
  using GetDocSLP_CompactLOUDS = dret::rmq::GetDocSLP<ExternalGenericStorage,
                                                       dret::Alphabet<>::int_width,
                                                       GCDASLP_CompactLOUDS>;
  using GetDocSLP_CSLP         = dret::rmq::GetDocSLP<ExternalGenericStorage,
                                                       dret::Alphabet<>::int_width,
                                                       GCDASLP_CSLP>;

  using SadaIdx = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                           dret::Alphabet<>,
                                           TCountIdx,
                                           dret::rmq::SadaCore<ExternalGenericStorage>>;
  using IlcpIdx = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                           dret::Alphabet<>,
                                           TCountIdx,
                                           dret::rmq::IlcpCore<ExternalGenericStorage>>;
  using CilcpIdx = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                            dret::Alphabet<>, TCountIdx, dret::rmq::CilcpCore<ExternalGenericStorage>>;

  // The GetDocSLP TSLP default intentionally matches gcda::DocListIdxGCDA<>::TSLP;
  // cache sharing depends on that exact type match because SDSL adds a type hash.
  using SadaIdxSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                              dret::Alphabet<>,
                                              TCountIdx,
                                              dret::rmq::SadaCore<ExternalGenericStorage,
                                                                  dret::Alphabet<>::int_width,
                                                                  sdsl::rmq_succinct_sct<true>,
                                                                  sdsl::sd_vector<>,
                                                                  GetDocSLP>>;
  using IlcpIdxSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                              dret::Alphabet<>,
                                              TCountIdx,
                                              dret::rmq::IlcpCore<ExternalGenericStorage,
                                                                  dret::Alphabet<>::int_width,
                                                                  sdsl::sd_vector<>,
                                                                  sdsl::rmq_succinct_sct<true>,
                                                                  sdsl::sd_vector<>,
                                                                  GetDocSLP>>;
  using CilcpIdxSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                               dret::Alphabet<>,
                                               TCountIdx,
                                               dret::rmq::CilcpCore<ExternalGenericStorage,
                                                                    dret::Alphabet<>::int_width,
                                                                    sdsl::sd_vector<>,
                                                                    sdsl::rmq_succinct_sct<true>,
                                                                    sdsl::sd_vector<>,
                                                                    GetDocSLP>>;

  // SADA / ILCP / CILCP-SLP variants for the compact GCDA TSLPs. Each consumes
  // the corresponding GCDA build's SLP cache (key match enforced by the
  // matching `TSLP` in `GetDocSLP_*`).
  template <typename TGetDoc>
  using SadaIdxSLPVariant = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                                     dret::Alphabet<>,
                                                     TCountIdx,
                                                     dret::rmq::SadaCore<ExternalGenericStorage,
                                                                         dret::Alphabet<>::int_width,
                                                                         sdsl::rmq_succinct_sct<true>,
                                                                         sdsl::sd_vector<>,
                                                                         TGetDoc>>;
  template <typename TGetDoc>
  using IlcpIdxSLPVariant = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                                     dret::Alphabet<>,
                                                     TCountIdx,
                                                     dret::rmq::IlcpCore<ExternalGenericStorage,
                                                                         dret::Alphabet<>::int_width,
                                                                         sdsl::sd_vector<>,
                                                                         sdsl::rmq_succinct_sct<true>,
                                                                         sdsl::sd_vector<>,
                                                                         TGetDoc>>;
  template <typename TGetDoc>
  using CilcpIdxSLPVariant = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                                      dret::Alphabet<>,
                                                      TCountIdx,
                                                      dret::rmq::CilcpCore<ExternalGenericStorage,
                                                                           dret::Alphabet<>::int_width,
                                                                           sdsl::sd_vector<>,
                                                                           sdsl::rmq_succinct_sct<true>,
                                                                           sdsl::sd_vector<>,
                                                                           TGetDoc>>;

  // SADA / ILCP / CILCP-SLP-NS: RMQ listing cores backed by the bare grammar::SLP<>
  // cache (kSLPNS) shared with dret::DocListIdxSLP. No bs/sf knobs.
  using SadaIdxSLP_NS = SadaIdxSLPVariant<GetDocSLP_NS>;
  using IlcpIdxSLP_NS = IlcpIdxSLPVariant<GetDocSLP_NS>;
  using CilcpIdxSLP_NS = CilcpIdxSLPVariant<GetDocSLP_NS>;

  // Bare-SLP container variants. Each typed grammar::SLP gets its own SDSL
  // type-hash, so cache files don't collide. GetDocSLP_NS_{Raw,DV,VV} share
  // their cache file with the corresponding DocListIdxSLP-{Raw,DV,VV} via
  // type matching.
  using BareSLP_Raw = grammar::SLP<>;
  using BareSLP_DV = grammar::SLP<sdsl::dac_vector<>, sdsl::dac_vector<>>;
  using BareSLP_VV = grammar::SLP<sdsl::vlc_vector<>, sdsl::vlc_vector<>>;
  using GetDocSLP_NS_Raw = dret::rmq::GetDocSLP_NS<ExternalGenericStorage,
                                                    dret::Alphabet<>::int_width,
                                                    BareSLP_Raw>;
  using GetDocSLP_NS_DV = dret::rmq::GetDocSLP_NS<ExternalGenericStorage,
                                                   dret::Alphabet<>::int_width,
                                                   BareSLP_DV>;
  using GetDocSLP_NS_VV = dret::rmq::GetDocSLP_NS<ExternalGenericStorage,
                                                   dret::Alphabet<>::int_width,
                                                   BareSLP_VV>;
  using SadaIdxSLP_NS_Raw  = SadaIdxSLPVariant<GetDocSLP_NS_Raw>;
  using IlcpIdxSLP_NS_Raw  = IlcpIdxSLPVariant<GetDocSLP_NS_Raw>;
  using CilcpIdxSLP_NS_Raw = CilcpIdxSLPVariant<GetDocSLP_NS_Raw>;
  using SadaIdxSLP_NS_DV  = SadaIdxSLPVariant<GetDocSLP_NS_DV>;
  using IlcpIdxSLP_NS_DV  = IlcpIdxSLPVariant<GetDocSLP_NS_DV>;
  using CilcpIdxSLP_NS_DV = CilcpIdxSLPVariant<GetDocSLP_NS_DV>;
  using SadaIdxSLP_NS_VV  = SadaIdxSLPVariant<GetDocSLP_NS_VV>;
  using IlcpIdxSLP_NS_VV  = IlcpIdxSLPVariant<GetDocSLP_NS_VV>;
  using CilcpIdxSLP_NS_VV = CilcpIdxSLPVariant<GetDocSLP_NS_VV>;

  // The GetDocDSLP default intentionally matches dgcda::DocListIdxDGCDA<>::TSLP.
  using SadaIdxDSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                               dret::Alphabet<>,
                                               TCountIdx,
                                               dret::rmq::SadaCore<ExternalGenericStorage,
                                                                   dret::Alphabet<>::int_width,
                                                                   sdsl::rmq_succinct_sct<true>,
                                                                   sdsl::sd_vector<>,
                                                                   GetDocDSLP>>;
  using IlcpIdxDSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                               dret::Alphabet<>,
                                               TCountIdx,
                                               dret::rmq::IlcpCore<ExternalGenericStorage,
                                                                   dret::Alphabet<>::int_width,
                                                                   sdsl::sd_vector<>,
                                                                   sdsl::rmq_succinct_sct<true>,
                                                                   sdsl::sd_vector<>,
                                                                   GetDocDSLP>>;
  using CilcpIdxDSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                                dret::Alphabet<>,
                                                TCountIdx,
                                                dret::rmq::CilcpCore<ExternalGenericStorage,
                                                                     dret::Alphabet<>::int_width,
                                                                     sdsl::sd_vector<>,
                                                                     sdsl::rmq_succinct_sct<true>,
                                                                     sdsl::sd_vector<>,
                                                                     GetDocDSLP>>;

  struct Config {
    IndexEnum index_t;
    std::size_t sampling_size = 0;
    uint32_t block_size = 512;
    float storing_factor = 4;
    GetDocEnum get_doc = GetDocEnum::DA;
    GCDASLPVariant gcda_slp = GCDASLPVariant::Default;
    BareSLPVariant bare_slp = BareSLPVariant::Default;
    PDLVariant pdl_variant = PDLVariant::Plain;
    PDLStoragePolicy pdl_storage_policy = PDLStoragePolicy::OccurrenceWeighted;

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
      return pdl_storage_policy < t_c.pdl_storage_policy;
    }
  };

  explicit Factory(dret::Config t_config, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : config_{std::move(t_config)} {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_));
    seq_size_ = buf.size();

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
    return MakeInner(t_config);
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
        // construct() would try to rebuild them from scratch (signature
        // mismatch on sr-index sampling).
        auto idx = std::make_shared<
            dret::DocListIdxBrute<ExternalGenericStorage, dret::Alphabet<>, sri::RIndex<ExternalGenericStorage>>>(
            std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::BRUTE_SR_INDEX: {
        auto idx = std::make_shared<dret::DocListIdxBrute<ExternalGenericStorage,
                                                          dret::Alphabet<>,
                                                          sri::SrIndexValidArea<ExternalGenericStorage>>>(
            std::ref(storage_), sri::SrIndexValidArea<ExternalGenericStorage>(storage_, t_config.sampling_size));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::GCDA: {
        switch (t_config.gcda_slp) {
          case GCDASLPVariant::CompactBP: {
            auto idx = std::make_shared<GCDAVariant<GCDASLP_CompactBP>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case GCDASLPVariant::CompactLOUDS: {
            auto idx = std::make_shared<GCDAVariant<GCDASLP_CompactLOUDS>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case GCDASLPVariant::CSLP: {
            auto idx = std::make_shared<GCDAVariant<GCDASLP_CSLP>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case GCDASLPVariant::Default:
          default: {
            auto idx = std::make_shared<dret::gcda::DocListIdxGCDA<ExternalGenericStorage>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
        }
        break;
      }

      case IndexEnum::DGCDA: {
        auto idx = std::make_shared<dret::dgcda::DocListIdxDGCDA<ExternalGenericStorage>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_OTF: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_OTF>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_CRL: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_CRL>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_EV: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_EV>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_DV: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_DV>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_VV: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_VV>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SADA: {
        if (t_config.get_doc == GetDocEnum::SLP_NS) {
          switch (t_config.bare_slp) {
            case BareSLPVariant::Raw: {
              typename SadaIdxSLP_NS_Raw::Core core(std::ref(storage_));
              auto idx = std::make_shared<SadaIdxSLP_NS_Raw>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::DV: {
              typename SadaIdxSLP_NS_DV::Core core(std::ref(storage_));
              auto idx = std::make_shared<SadaIdxSLP_NS_DV>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::VV: {
              typename SadaIdxSLP_NS_VV::Core core(std::ref(storage_));
              auto idx = std::make_shared<SadaIdxSLP_NS_VV>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::Default:
            default: {
              typename SadaIdxSLP_NS::Core core(std::ref(storage_));
              auto idx = std::make_shared<SadaIdxSLP_NS>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
          }
          break;
        }
        if (t_config.get_doc == GetDocEnum::SLP) {
          switch (t_config.gcda_slp) {
            case GCDASLPVariant::CompactBP: {
              using IdxT = SadaIdxSLPVariant<GetDocSLP_CompactBP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CompactLOUDS: {
              using IdxT = SadaIdxSLPVariant<GetDocSLP_CompactLOUDS>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CSLP: {
              using IdxT = SadaIdxSLPVariant<GetDocSLP_CSLP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::Default:
            default: {
              typename SadaIdxSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<SadaIdxSLP>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
          }
          break;
        }
        if (t_config.get_doc == GetDocEnum::DSLP) {
          typename SadaIdxDSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
          auto idx = std::make_shared<SadaIdxDSLP>(std::ref(storage_), core);
          construct(*idx, config_);
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
          break;
        }

        auto idx = std::make_shared<SadaIdx>(std::ref(storage_));
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::ILCP: {
        if (t_config.get_doc == GetDocEnum::SLP_NS) {
          switch (t_config.bare_slp) {
            case BareSLPVariant::Raw: {
              typename IlcpIdxSLP_NS_Raw::Core core(std::ref(storage_));
              auto idx = std::make_shared<IlcpIdxSLP_NS_Raw>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::DV: {
              typename IlcpIdxSLP_NS_DV::Core core(std::ref(storage_));
              auto idx = std::make_shared<IlcpIdxSLP_NS_DV>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::VV: {
              typename IlcpIdxSLP_NS_VV::Core core(std::ref(storage_));
              auto idx = std::make_shared<IlcpIdxSLP_NS_VV>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::Default:
            default: {
              typename IlcpIdxSLP_NS::Core core(std::ref(storage_));
              auto idx = std::make_shared<IlcpIdxSLP_NS>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
          }
          break;
        }
        if (t_config.get_doc == GetDocEnum::SLP) {
          switch (t_config.gcda_slp) {
            case GCDASLPVariant::CompactBP: {
              using IdxT = IlcpIdxSLPVariant<GetDocSLP_CompactBP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CompactLOUDS: {
              using IdxT = IlcpIdxSLPVariant<GetDocSLP_CompactLOUDS>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CSLP: {
              using IdxT = IlcpIdxSLPVariant<GetDocSLP_CSLP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::Default:
            default: {
              typename IlcpIdxSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IlcpIdxSLP>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
          }
          break;
        }
        if (t_config.get_doc == GetDocEnum::DSLP) {
          typename IlcpIdxDSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
          auto idx = std::make_shared<IlcpIdxDSLP>(std::ref(storage_), core);
          construct(*idx, config_);
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
          break;
        }

        auto idx = std::make_shared<IlcpIdx>(std::ref(storage_));
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::PDL: {
        auto policy = toPDLStoragePolicy(t_config.pdl_storage_policy);
        // Lambda: given the concrete TIndex, instantiate & load.
        auto build = [this, &t_config, policy, &index](auto type_tag) {
          using TIndex = typename decltype(type_tag)::type;
          auto idx = std::make_shared<TIndex>(
              std::ref(storage_),
              t_config.block_size,
              t_config.storing_factor,
              policy);
          construct(*idx, config_);
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
        };
        // Tag helper so the lambda can deduce a concrete type.
        struct Plain_DA   { using type = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage>; };
        struct Plain_SLP  { using type = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage,
                                                                       dret::Alphabet<>,
                                                                       TCountIdx,
                                                                       dret::pdl::PDLGetDocsSLP<ExternalGenericStorage>>; };
        struct Plain_DSLP { using type = dret::pdl::DocListIdxPDLPlain<ExternalGenericStorage,
                                                                       dret::Alphabet<>,
                                                                       TCountIdx,
                                                                       dret::pdl::PDLGetDocsDSLP<ExternalGenericStorage>>; };
        struct RP_DA      { using type = dret::pdl::DocListIdxPDLRP<ExternalGenericStorage>; };
        struct RP_SLP     { using type = dret::pdl::DocListIdxPDLRP<ExternalGenericStorage,
                                                                    dret::Alphabet<>,
                                                                    TCountIdx,
                                                                    dret::pdl::PDLGetDocsSLP<ExternalGenericStorage>>; };
        struct RP_DSLP    { using type = dret::pdl::DocListIdxPDLRP<ExternalGenericStorage,
                                                                    dret::Alphabet<>,
                                                                    TCountIdx,
                                                                    dret::pdl::PDLGetDocsDSLP<ExternalGenericStorage>>; };
        struct BC_DA      { using type = dret::pdl::DocListIdxPDLBC<ExternalGenericStorage>; };
        struct BC_SLP     { using type = dret::pdl::DocListIdxPDLBC<ExternalGenericStorage,
                                                                    dret::Alphabet<>,
                                                                    TCountIdx,
                                                                    dret::pdl::PDLGetDocsSLP<ExternalGenericStorage>>; };
        struct BC_DSLP    { using type = dret::pdl::DocListIdxPDLBC<ExternalGenericStorage,
                                                                    dret::Alphabet<>,
                                                                    TCountIdx,
                                                                    dret::pdl::PDLGetDocsDSLP<ExternalGenericStorage>>; };

        switch (t_config.pdl_variant) {
          case PDLVariant::Plain:
            switch (t_config.get_doc) {
              case GetDocEnum::SLP:  build(Plain_SLP{});  break;
              case GetDocEnum::DSLP: build(Plain_DSLP{}); break;
              case GetDocEnum::DA:
              default:               build(Plain_DA{});   break;
            }
            break;
          case PDLVariant::RP:
            switch (t_config.get_doc) {
              case GetDocEnum::SLP:  build(RP_SLP{});  break;
              case GetDocEnum::DSLP: build(RP_DSLP{}); break;
              case GetDocEnum::DA:
              default:               build(RP_DA{});   break;
            }
            break;
          case PDLVariant::BC:
            switch (t_config.get_doc) {
              case GetDocEnum::SLP:  build(BC_SLP{});  break;
              case GetDocEnum::DSLP: build(BC_DSLP{}); break;
              case GetDocEnum::DA:
              default:               build(BC_DA{});   break;
            }
            break;
        }
        break;
      }

      case IndexEnum::SLP_NS: {
        switch (t_config.bare_slp) {
          case BareSLPVariant::Raw: {
            auto idx = std::make_shared<dret::DocListIdxSLP<ExternalGenericStorage,
                                                            dret::Alphabet<>,
                                                            sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                                                            BareSLP_Raw>>(std::ref(storage_));
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case BareSLPVariant::DV: {
            auto idx = std::make_shared<dret::DocListIdxSLP<ExternalGenericStorage,
                                                            dret::Alphabet<>,
                                                            sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                                                            BareSLP_DV>>(std::ref(storage_));
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case BareSLPVariant::VV: {
            auto idx = std::make_shared<dret::DocListIdxSLP<ExternalGenericStorage,
                                                            dret::Alphabet<>,
                                                            sri::RIndexCount<ExternalGenericStorage, dret::Alphabet<>>,
                                                            BareSLP_VV>>(std::ref(storage_));
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case BareSLPVariant::Default:
          default: {
            auto idx = std::make_shared<dret::DocListIdxSLP<ExternalGenericStorage>>(std::ref(storage_));
            construct(*idx, config_);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
        }
        break;
      }

      case IndexEnum::CILCP: {
        if (t_config.get_doc == GetDocEnum::SLP_NS) {
          switch (t_config.bare_slp) {
            case BareSLPVariant::Raw: {
              typename CilcpIdxSLP_NS_Raw::Core core(std::ref(storage_));
              auto idx = std::make_shared<CilcpIdxSLP_NS_Raw>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::DV: {
              typename CilcpIdxSLP_NS_DV::Core core(std::ref(storage_));
              auto idx = std::make_shared<CilcpIdxSLP_NS_DV>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::VV: {
              typename CilcpIdxSLP_NS_VV::Core core(std::ref(storage_));
              auto idx = std::make_shared<CilcpIdxSLP_NS_VV>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case BareSLPVariant::Default:
            default: {
              typename CilcpIdxSLP_NS::Core core(std::ref(storage_));
              auto idx = std::make_shared<CilcpIdxSLP_NS>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
          }
          break;
        }
        if (t_config.get_doc == GetDocEnum::SLP) {
          switch (t_config.gcda_slp) {
            case GCDASLPVariant::CompactBP: {
              using IdxT = CilcpIdxSLPVariant<GetDocSLP_CompactBP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CompactLOUDS: {
              using IdxT = CilcpIdxSLPVariant<GetDocSLP_CompactLOUDS>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CSLP: {
              using IdxT = CilcpIdxSLPVariant<GetDocSLP_CSLP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::Default:
            default: {
              typename CilcpIdxSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<CilcpIdxSLP>(std::ref(storage_), core);
              construct(*idx, config_);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
          }
          break;
        }
        if (t_config.get_doc == GetDocEnum::DSLP) {
          typename CilcpIdxDSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
          auto idx = std::make_shared<CilcpIdxDSLP>(std::ref(storage_), core);
          construct(*idx, config_);
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
          break;
        }

        auto idx = std::make_shared<CilcpIdx>(std::ref(storage_));
        construct(*idx, config_);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
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
