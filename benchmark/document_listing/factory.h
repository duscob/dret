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
#include "sr-index/sr_idx_generic.h"
#include "sr-index/sr_index.h"

#include "../tool/definitions.h"

#include "dret/basic_slp_span_length.h"
#include "dret/config.h"
#include "dret/differential_light_slp.h"
#include "dret/doc_list_idx_slp.h"
#include "dret/doc_list_index.h"
#include "dret/doc_list_index_brute.h"
#include "dret/doc_list_index_rmq.h"
#include "dret/doc_list_sampled_tree_dgcda.h"
#include "dret/doc_list_sampled_tree_gcda.h"


using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

template <uint8_t t_width = 8>
class Factory {
 public:
  enum class IndexEnum {
    BRUTE_R_INDEX,
    BRUTE_SR_INDEX,
    GCDA,
    DGCDA,
    DGCDA_OTF,  // BasicSLPOnTheFlySpanLength — no stored lengths
    DGCDA_CRL,  // BasicSLPCachedRootSpanLengths — lengths cached for roots only
    DGCDA_EV,   // default SLP; TRoots/TSpanSums/TSamples = sdsl::enc_vector<>
    DGCDA_DV,   // default SLP; TRoots/TSpanSums/TSamples = sdsl::dac_vector<>
    DGCDA_VV,   // default SLP; TRoots/TSpanSums/TSamples = sdsl::vlc_vector<>
    SADA,       // RMinQ on prev_doc
    ILCP,       // RMinQ on backward-ILCP runs
    CILCP,      // RMinQ on doc-aware compressed backward-ILCP runs
    SLP_NS,     // Phase C: dret::DocListIdxSLP — non-sampled grammar::SLP<>
  };

  enum class GetDocEnum {
    DA,
    SLP,
    DSLP,
  };

  // GCDA's TSLP choice — independent from the RMQ DA-lookup `GetDocEnum`. The
  // RMQ-SLP path consults this too so SADA/ILCP/CILCP-SLP reuse the SLP cache
  // file built for the matching GCDA variant.
  enum class GCDASLPVariant {
    Default,       // grammar::LightSLP<...> — the existing GCDA default
    CompactBP,     // grammar::CompactBPSLP<>
    CompactLOUDS,  // grammar::CompactLOUDSSLP<>
    CSLP,          // grammar::CombinedSLPWithUnitCover<> — Phase B
  };

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
      sri::SrIdxGeneric<sri::SrIndexValidArea<ExternalGenericStorage, dret::Alphabet<>>, 16>,
      TSLP>;

  // GCDA TSLP variants (Phase A — Compact-BP / Compact-LOUDS; Phase B — CSLP).
  using GCDASLP_CompactBP    = grammar::CompactBPSLP<>;
  using GCDASLP_CompactLOUDS = grammar::CompactLOUDSSLP<>;
  using GCDASLP_CSLP         = grammar::CombinedSLPWithUnitCover<>;

  template <typename TSLP>
  using GCDAVariant = dret::gcda::DocListIdxGCDA<
      ExternalGenericStorage,
      dret::Alphabet<>,
      sri::SrIdxGeneric<sri::SrIndexValidArea<ExternalGenericStorage, dret::Alphabet<>>, 16>,
      TSLP>;

  using TCountIdx = sri::SrIdxGeneric<sri::SrIndexValidArea<ExternalGenericStorage, dret::Alphabet<>>, 16>;
  using GetDocSLP = dret::rmq::GetDocSLP<ExternalGenericStorage>;
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
                                                                  sdsl::bit_vector,
                                                                  sdsl::rmq_succinct_sct<true>,
                                                                  sdsl::sd_vector<>,
                                                                  GetDocSLP>>;
  using CilcpIdxSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                               dret::Alphabet<>,
                                               TCountIdx,
                                               dret::rmq::CilcpCore<ExternalGenericStorage,
                                                                    dret::Alphabet<>::int_width,
                                                                    sdsl::bit_vector,
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
                                                                         sdsl::bit_vector,
                                                                         sdsl::rmq_succinct_sct<true>,
                                                                         sdsl::sd_vector<>,
                                                                         TGetDoc>>;
  template <typename TGetDoc>
  using CilcpIdxSLPVariant = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                                      dret::Alphabet<>,
                                                      TCountIdx,
                                                      dret::rmq::CilcpCore<ExternalGenericStorage,
                                                                           dret::Alphabet<>::int_width,
                                                                           sdsl::bit_vector,
                                                                           sdsl::rmq_succinct_sct<true>,
                                                                           sdsl::sd_vector<>,
                                                                           TGetDoc>>;

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
                                                                   sdsl::bit_vector,
                                                                   sdsl::rmq_succinct_sct<true>,
                                                                   sdsl::sd_vector<>,
                                                                   GetDocDSLP>>;
  using CilcpIdxDSLP = dret::rmq::DocListIdxRMQ<ExternalGenericStorage,
                                                dret::Alphabet<>,
                                                TCountIdx,
                                                dret::rmq::CilcpCore<ExternalGenericStorage,
                                                                     dret::Alphabet<>::int_width,
                                                                     sdsl::bit_vector,
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
      return gcda_slp < t_c.gcda_slp;
    }
  };

  explicit Factory(dret::Config t_config, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : config_{std::move(t_config)} {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_));
    seq_size_ = buf.size();

    r_index_ = std::make_shared<sri::RIndex<ExternalGenericStorage>>(std::ref(storage_));
    r_index_->load(config_);

    sr_index_ = std::make_shared<sri::SrIndexValidArea<ExternalGenericStorage>>(std::ref(storage_), 8);
    sr_index_->load(config_);

    load(doc_endings_);
    load(doc_endings_rank_, [this]() {
      return TDocEndingsRank(&this->doc_endings_.item);
    });
    n_doc_ = doc_endings_rank_.item(doc_endings_.item.size());
  }

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

  Index MakeIndex(const Config& t_config) {
    auto it = indexes_.find(t_config);
    if (it != indexes_.end()) {
      return it->second;
    }

    Index index;
    switch (t_config.index_t) {
      case IndexEnum::BRUTE_R_INDEX: {
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
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case GCDASLPVariant::CompactLOUDS: {
            auto idx = std::make_shared<GCDAVariant<GCDASLP_CompactLOUDS>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case GCDASLPVariant::CSLP: {
            auto idx = std::make_shared<GCDAVariant<GCDASLP_CSLP>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
            idx->load(config_);
            index = {idx, sdsl::size_in_bytes(*idx)};
            break;
          }
          case GCDASLPVariant::Default:
          default: {
            auto idx = std::make_shared<dret::gcda::DocListIdxGCDA<ExternalGenericStorage>>(
                std::ref(storage_), t_config.block_size, t_config.storing_factor);
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
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_OTF: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_OTF>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_CRL: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_CRL>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_EV: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_EV>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_DV: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_DV>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::DGCDA_VV: {
        auto idx = std::make_shared<DGCDAVariant<DGCDASLP_VV>>(
            std::ref(storage_), t_config.block_size, t_config.storing_factor);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SADA: {
        if (t_config.get_doc == GetDocEnum::SLP) {
          switch (t_config.gcda_slp) {
            case GCDASLPVariant::CompactBP: {
              using IdxT = SadaIdxSLPVariant<GetDocSLP_CompactBP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CompactLOUDS: {
              using IdxT = SadaIdxSLPVariant<GetDocSLP_CompactLOUDS>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CSLP: {
              using IdxT = SadaIdxSLPVariant<GetDocSLP_CSLP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::Default:
            default: {
              typename SadaIdxSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<SadaIdxSLP>(std::ref(storage_), core);
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
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
          break;
        }

        auto idx = std::make_shared<SadaIdx>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::ILCP: {
        if (t_config.get_doc == GetDocEnum::SLP) {
          switch (t_config.gcda_slp) {
            case GCDASLPVariant::CompactBP: {
              using IdxT = IlcpIdxSLPVariant<GetDocSLP_CompactBP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CompactLOUDS: {
              using IdxT = IlcpIdxSLPVariant<GetDocSLP_CompactLOUDS>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CSLP: {
              using IdxT = IlcpIdxSLPVariant<GetDocSLP_CSLP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::Default:
            default: {
              typename IlcpIdxSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IlcpIdxSLP>(std::ref(storage_), core);
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
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
          break;
        }

        auto idx = std::make_shared<IlcpIdx>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SLP_NS: {
        auto idx = std::make_shared<dret::DocListIdxSLP<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::CILCP: {
        if (t_config.get_doc == GetDocEnum::SLP) {
          switch (t_config.gcda_slp) {
            case GCDASLPVariant::CompactBP: {
              using IdxT = CilcpIdxSLPVariant<GetDocSLP_CompactBP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CompactLOUDS: {
              using IdxT = CilcpIdxSLPVariant<GetDocSLP_CompactLOUDS>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::CSLP: {
              using IdxT = CilcpIdxSLPVariant<GetDocSLP_CSLP>;
              typename IdxT::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<IdxT>(std::ref(storage_), core);
              idx->load(config_);
              index = {idx, sdsl::size_in_bytes(*idx)};
              break;
            }
            case GCDASLPVariant::Default:
            default: {
              typename CilcpIdxSLP::Core core(std::ref(storage_), t_config.block_size, t_config.storing_factor);
              auto idx = std::make_shared<CilcpIdxSLP>(std::ref(storage_), core);
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
          idx->load(config_);
          index = {idx, sdsl::size_in_bytes(*idx)};
          break;
        }

        auto idx = std::make_shared<CilcpIdx>(std::ref(storage_));
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
        auto locate = [this](const auto& tt_pattern) {
          return this->r_index_->Locate(tt_pattern);
        };

        index = new dret::DocListIndexBrute(locate, doc_endings_rank_.item);
        index_size = sdsl::size_in_bytes(*r_index_) + doc_endings_rank_.size_in_bytes;
        break;
      }

      case IndexEnum::BRUTE_SR_INDEX: {
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
