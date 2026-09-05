//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/9/25.
//

#pragma once

#include "dret/repair.h"

#include "sr-index/alphabet.h"
#include "sr-index/config.h"
#include "sr-index/index_base.h"

namespace dret {

using sri::GenericStorage;
using sri::get;
using sri::set;

using sri::JSON;

using sri::Alphabet;

//~~~~~~~

namespace conf {
using namespace sri::conf;
constexpr std::string_view kDocEnds = "docEnd";
constexpr std::string_view kDA = "da";
constexpr std::string_view kGCDA = "gcda";
constexpr std::string_view kSLP = "slp";
constexpr std::string_view kDocs = "docs";
// The SLP over the document array plus its compact sequence, as parsed from
// irepair's .R/.C output. Both depend only on the DA, so unlike kSLP/kDocs
// these are cached WITHOUT the "<bs>-<sf>_" prefix and are shared by every cell
// of a (block size, storing factor) sweep. Before they existed GCDA re-parsed
// the grammar in all 20 cells of the grid; DGCDA already did the equivalent via
// its own dgcda_slp_grammar.
constexpr std::string_view kSLPGrammar = "slpGrammar";
constexpr std::string_view kSLPCompactSeq = "slpCompactSeq";
constexpr std::string_view kDGCDA = "dgcda";
// Phase C: non-sampled grammar::SLP<> cache (dret::DocListIdxSLP).
// Distinct from kGCDA::kSLP so the bare SLP and the GCDA-sampled SLPs
// don't share a logical key (their type-hashes are already distinct,
// but a separate prefix makes the on-disk files easier to inspect).
constexpr std::string_view kSLPNS = "slpNS";
// Non-sampled *differential* SLP cache (bare-diff): base dret::DifferentialSLP,
// no sampled tree / GCChunks. Distinct key from the (bs/sf-keyed) DGCDA
// DifferentialLightSLP cache — bare-diff has a single fixed internal sample
// block_size and no storing_factor.
constexpr std::string_view kDSLPNS = "dslpNS";

constexpr std::string_view kSADA = "sada";
constexpr std::string_view kILCP = "ilcp";
constexpr std::string_view kCILCP = "cilcp";
constexpr std::string_view kRmq = "rmq";
constexpr std::string_view kRunHeads = "runHeads";
constexpr std::string_view kRunValues = "runValues";
constexpr std::string_view kPrevDoc = "prevDoc";
constexpr std::string_view kRmqNDoc = "rmq_n_doc";
// Backward interleaved-LCP array (ComputeIlcpBackward). Collection-level and
// shared by every ILCP-family core (ILCP / CILCP and their -L variants), so it
// is cached once instead of recomputed per core.
constexpr std::string_view kIlcpArray = "ilcp_array";

// PDL (precomputed document listing) — sparse suffix-tree indexes with
// per-variant stored-set codecs. kPDL is the umbrella; kTree is the
// shared tree-topology cache; each variant (kPlain/kRP/kBC) owns its
// own kSets payload, and kBC additionally owns kDict.
// RLCSA sidecar cache (Track C of pdl_rlcsa_baseline_plan.md).
// RLCSA writes its own on-disk files via writeTo(base); the key maps to the
// file basename, not an SDSL typed-cache entry.
constexpr std::string_view kRLCSA = "rlcsa";

constexpr std::string_view kPDL = "pdl";
constexpr std::string_view kPlain = "plain";
constexpr std::string_view kRP = "rp";
constexpr std::string_view kBC = "bc";
constexpr std::string_view kTree = "tree";
constexpr std::string_view kSets = "sets";
constexpr std::string_view kDict = "dict";
}  // namespace conf

template <uint8_t t_width = DRET_DEFAULT_ALPHABET_WIDTH>
struct Keys {
  Keys() {
    keys = sri::createDefaultKeys<t_width>();

    keys.update({
        {conf::kDocEnds, "doc_end"},
        {conf::kDA, "da"},
        {conf::kSLPNS, "slp_ns"},
        {conf::kDSLPNS, "dslp_ns"},
        {
            conf::kGCDA,
            {
                {conf::kSLP, "gcda_slp"},
                {conf::kDocs, "gcda_docs"},
                {conf::kSLPGrammar, "gcda_slp_grammar"},
                {conf::kSLPCompactSeq, "gcda_slp_compact_seq"},
            },
        },
        {
            conf::kDGCDA,
            {
                {conf::kSLP, "dgcda_slp"},
                {conf::kDocs, "dgcda_docs"},
            },
        },
        // One namespace per RMQ family, holding every artefact that family can
        // need. Each core takes the subset it uses: the -L cores stop at the
        // partition (rmq, run_heads), the published cores additionally read the
        // value array (prev_doc / run_values) that their value-based stop
        // consults. There is no separate namespace for the published cores --
        // they are the base case, not a variant of it.
        {
            conf::kSADA,
            {
                {conf::kRmq, "sada_rmq"},
                {conf::kPrevDoc, "sada_prev_doc"},
            },
        },
        {
            conf::kILCP,
            {
                {conf::kRmq, "ilcp_rmq"},
                {conf::kRunHeads, "ilcp_run_heads"},
                {conf::kRunValues, "ilcp_run_values"},
            },
        },
        {
            conf::kCILCP,
            {
                {conf::kRmq, "cilcp_rmq"},
                {conf::kRunHeads, "cilcp_run_heads"},
                {conf::kRunValues, "cilcp_run_values"},
            },
        },
        {conf::kRLCSA, "rlcsa"},
        {conf::kRmqNDoc, "rmq_n_doc"},
        {conf::kIlcpArray, "ilcp_array"},
        {
            conf::kPDL,
            {
                {conf::kTree, "pdl_tree"},
                {
                    conf::kPlain,
                    {
                        {conf::kSets, "pdl_plain_sets"},
                    },
                },
                {
                    conf::kRP,
                    {
                        {conf::kSets, "pdl_rp_sets"},
                    },
                },
                {
                    conf::kBC,
                    {
                        {conf::kSets, "pdl_bc_sets"},
                        {conf::kDict, "pdl_bc_dict"},
                    },
                },
            },
        },
    });
  }

  nlohmann::json keys;
};

const Keys<> kDefaultKeys;

template <uint8_t t_width>
const auto& createDefaultKeys() {
  static Keys<t_width> keys;
  return keys.keys;
}

struct Config : public sri::Config {
  uint64_t data_delim = 0;

  // How to run irepair over this collection's document array. Defaults derive
  // everything from the array's size, so this needs setting only when a host
  // cannot afford the derived value. It lives here, rather than being read from
  // the environment, because <MB> changes the grammar that gets built: the value
  // used has to travel with the run's configuration.
  repair::Options repair;

  Config() = default;

  Config(const std::filesystem::path& t_data_path,
         const std::filesystem::path& t_output_dir,
         sri::SAAlgo t_sa_algo,
         bool t_delete_files = false,
         uint8_t t_data_width = 8,
         uint64_t t_data_delim = 0,
         JSON t_keys = kDefaultKeys.keys)
      : sri::Config(t_data_path, t_output_dir, t_sa_algo, t_delete_files, t_data_width, std::move(t_keys)),
        data_delim(t_data_delim) {}
};

}  // namespace dret
