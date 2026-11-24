//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/9/25.
//

#pragma once

#include <format>
#include <string>
#include <string_view>

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
// Cache components of the document-array grammars, shared by every index that
// contains them (see dret/cache_components.h). Collection-level ones carry no
// prefix; the others carry the parameters their construction reads, through one
// of the prefix patterns below. Each name says what the file holds; its type
// hash says how it is encoded.
constexpr std::string_view kDaGrammar = "daGrammar";              // RePair grammar of the DA
constexpr std::string_view kDaGrammarSeq = "daGrammarSeq";        // its top-level sequence
constexpr std::string_view kDaCnfGrammar = "daCnfGrammar";        // the DA grammar in CNF
constexpr std::string_view kDaCnfGrammarBP = "daCnfGrammarBP";    // CNF grammar as a BP tree
constexpr std::string_view kDaCnfGrammarLOUDS = "daCnfGrammarLOUDS";
constexpr std::string_view kDaSampledTree = "daSampledTree";      // the sampled tree of the CNF grammar
constexpr std::string_view kDaSampledLeaves = "daSampledLeaves";  // its leaves as CNF variables
constexpr std::string_view kDaSampledLeavesBP = "daSampledLeavesBP";
constexpr std::string_view kDaSampledLeavesLOUDS = "daSampledLeavesLOUDS";
constexpr std::string_view kDaLeafCovers = "daLeafCovers";        // its leaves as covers in the DA grammar
constexpr std::string_view kDaNodeDocListsPlain = "daNodeDocListsPlain";  // document list per sampled node
constexpr std::string_view kDaNodeDocListsRP = "daNodeDocListsRP";
constexpr std::string_view kDaDiffGrammar = "daDiffGrammar";      // RePair grammar of the differential DA
constexpr std::string_view kDaDiffRoots = "daDiffRoots";          // its top-level sequence
constexpr std::string_view kDaDiffSpanSums = "daDiffSpanSums";    // differential sum of each rule
constexpr std::string_view kDaDiffSamples = "daDiffSamples";      // absolute samples every spacing positions

// PDL components: the collapsed suffix tree, the nodes the policy selects, and
// the document lists of those nodes in each codec.
constexpr std::string_view kPdlTree = "pdlTree";
constexpr std::string_view kPdlSelection = "pdlSelection";
constexpr std::string_view kPdlDocListsPlain = "pdlDocListsPlain";
constexpr std::string_view kPdlDocListsRP = "pdlDocListsRP";
constexpr std::string_view kPdlDocListsBC = "pdlDocListsBC";

// Prefix patterns for cache keys that carry parameters, each named after the
// parameters it carries. Every parameterised key is built from one of these
// through PrefixedKey, so the file names an index can produce are all declared
// here, next to the component names.
constexpr std::string_view kPrefix = "prefix";
constexpr std::string_view kBlkSf = "blkSf";          // block size, storing factor
constexpr std::string_view kSpc = "spc";              // sample spacing
constexpr std::string_view kBlk = "blk";              // block size
constexpr std::string_view kBlkPol = "blkPol";        // block size, selection policy
constexpr std::string_view kBlkSfPol = "blkSfPol";    // block size, storing factor, selection policy
// The term each PDL selection policy contributes to a key.
constexpr std::string_view kPolicy = "policy";
constexpr std::string_view kOccW = "occw";
constexpr std::string_view kListW = "listw";
constexpr std::string_view kLeaves = "leaves";
constexpr std::string_view kAllNodes = "all";

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

// RLCSA sidecar cache. RLCSA writes its own on-disk files via writeTo(base); the
// key maps to the file basename, not an SDSL typed-cache entry.
constexpr std::string_view kRLCSA = "rlcsa";

}  // namespace conf

template <uint8_t t_width = DRET_DEFAULT_ALPHABET_WIDTH>
struct Keys {
  Keys() {
    keys = sri::createDefaultKeys<t_width>();

    keys.update({
        {conf::kDocEnds, "doc_end"},
        {conf::kDA, "da"},
        {conf::kDaGrammar, "da_grammar"},
        {conf::kDaGrammarSeq, "da_grammar_seq"},
        {conf::kDaCnfGrammar, "da_cnf_grammar"},
        {conf::kDaCnfGrammarBP, "da_cnf_grammar_bp"},
        {conf::kDaCnfGrammarLOUDS, "da_cnf_grammar_louds"},
        {conf::kDaSampledTree, "da_sampled_tree"},
        {conf::kDaSampledLeaves, "da_sampled_leaves"},
        {conf::kDaSampledLeavesBP, "da_sampled_leaves_bp"},
        {conf::kDaSampledLeavesLOUDS, "da_sampled_leaves_louds"},
        {conf::kDaLeafCovers, "da_leaf_covers"},
        {conf::kDaNodeDocListsPlain, "da_node_doclists_plain"},
        {conf::kDaNodeDocListsRP, "da_node_doclists_rp"},
        {conf::kDaDiffGrammar, "da_diff_grammar"},
        {conf::kDaDiffRoots, "da_diff_roots"},
        {conf::kDaDiffSpanSums, "da_diff_span_sums"},
        {conf::kDaDiffSamples, "da_diff_samples"},
        {conf::kPdlTree, "pdl_tree"},
        {conf::kPdlSelection, "pdl_selection"},
        {conf::kPdlDocListsPlain, "pdl_doclists_plain"},
        {conf::kPdlDocListsRP, "pdl_doclists_rp"},
        {conf::kPdlDocListsBC, "pdl_doclists_bc"},
        {
            conf::kPrefix,
            {
                {conf::kBlkSf, "blk{}-sf{}_"},
                {conf::kSpc, "spc{}_"},
                {conf::kBlk, "blk{}_"},
                {conf::kBlkPol, "blk{}-{}_"},
                {conf::kBlkSfPol, "blk{}-sf{}-{}_"},
            },
        },
        {
            conf::kPolicy,
            {
                {conf::kOccW, "occw"},
                {conf::kListW, "listw"},
                {conf::kLeaves, "leaves"},
                {conf::kAllNodes, "all"},
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

// Storing factor infinity: no node's list is ever worth storing over its
// children's (GCDA) and no node's occurrences ever exceed it times its documents
// (PDL), so only the leaves keep lists. A finite sentinel, not a float infinity:
// the build uses -Ofast, whose fast-math assumes there are no infinities and folds
// any test for one. Every rule multiplies it by a list size, at most ~1e8, far
// below the float maximum (~3.4e38), so nothing overflows.
constexpr float kInfiniteStoringFactor = 1e30f;

constexpr bool IsInfiniteStoringFactor(float t_storing_factor) {
  return t_storing_factor >= kInfiniteStoringFactor;
}

// How a storing factor reads in cache keys and benchmark names: "inf", or the
// number as std::format prints it (so the keys of finite values are unchanged).
inline std::string StoringFactorTerm(float t_storing_factor) {
  return IsInfiniteStoringFactor(t_storing_factor) ? std::string("inf") : std::format("{}", t_storing_factor);
}

// Builds a parameterised cache key: the prefix pattern t_keys[prefix][t_scheme],
// formatted with t_args, followed by t_name. The only place such keys are built.
template <typename... TArgs>
std::string PrefixedKey(const JSON& t_keys, std::string_view t_scheme, const std::string& t_name,
                        const TArgs&... t_args) {
  const auto pattern = t_keys[conf::kPrefix][t_scheme].template get<std::string>();
  return std::vformat(pattern, std::make_format_args(t_args...)) + t_name;
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
