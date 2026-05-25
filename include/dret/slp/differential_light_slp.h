//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/22/26.
//

#pragma once

#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "dret/config.h"
#include "dret/slp/differential_slp.h"

namespace dret {

template <typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>,
          typename TSampledSLP = grammar::SampledSLP<>,
          typename TRoots = sdsl::int_vector<>,
          typename TSpanSums = sdsl::int_vector<>,
          typename TSamples = sdsl::int_vector<>,
          typename TSampleRootsPos = sdsl::enc_vector<>,
          typename TBV = sdsl::sd_vector<>>
class DifferentialLightSLP : public DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>,
                             public TSampledSLP {
  using DiffBase = DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>;

 public:
  using size_type = std::size_t;

  DifferentialLightSLP() = default;

  template <typename TOtherSLP,
            typename TOtherSampledSLP,
            typename TOtherRoots,
            typename TOtherSpanSums,
            typename TOtherSamples,
            typename TOtherSampleRootsPos,
            typename TOtherBV,
            typename ActionSLPRules = grammar::NoAction,
            typename ActionSLPLengths = grammar::NoAction,
            typename ActionIntContainers = grammar::NoAction>
  DifferentialLightSLP(const DifferentialLightSLP<TOtherSLP,
                                                  TOtherSampledSLP,
                                                  TOtherRoots,
                                                  TOtherSpanSums,
                                                  TOtherSamples,
                                                  TOtherSampleRootsPos,
                                                  TOtherBV>& other,
                       ActionSLPRules&& action_slp_rules = grammar::NoAction(),
                       ActionSLPLengths&& action_slp_lengths = grammar::NoAction(),
                       ActionIntContainers&& action_int_containers = grammar::NoAction())
      : DiffBase(static_cast<const DifferentialSLP<TOtherSLP,
                                                   TOtherRoots,
                                                   TOtherSpanSums,
                                                   TOtherSamples,
                                                   TOtherSampleRootsPos,
                                                   TOtherBV>&>(other),
                 std::forward<ActionSLPRules>(action_slp_rules),
                 std::forward<ActionSLPLengths>(action_slp_lengths),
                 std::forward<ActionIntContainers>(action_int_containers)),
        TSampledSLP(static_cast<const TOtherSampledSLP&>(other)) {}

  // cached_cnf: optional pre-built CNF SLP of the DA (RePair output) for the
  // sampled tree. It is bs/sf/variant-independent, so the construct path builds
  // it once per collection and reuses it across all cells, avoiding a redundant
  // RePair per (bs,sf). Pass nullptr to build it inline.
  void Compute(const sdsl::int_vector<>& da,
               uint32_t block_size,
               float storing_factor,
               grammar::Chunks<>& cslp_docs_out,
               const grammar::SLP<>* cached_cnf = nullptr) {
    DiffBase::Compute(da, block_size);
    buildSampledSLP(da, block_size, storing_factor, cslp_docs_out, cached_cnf);
  }

  // Build only the sampled-tree part (the DiffBase grammar/fields must already be
  // set, e.g. via DiffBase::RestoreBaseGrammar + FinishCompute). Lets the construct
  // path drive grammar caching for both the diff grammar and the sampled-tree CNF.
  void BuildSampled(const sdsl::int_vector<>& da,
                    uint32_t block_size,
                    float storing_factor,
                    grammar::Chunks<>& cslp_docs_out,
                    const grammar::SLP<>* cached_cnf) {
    buildSampledSLP(da, block_size, storing_factor, cslp_docs_out, cached_cnf);
  }

  std::size_t serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, const std::string& name = "") const {
    return DiffBase::serialize(out) + TSampledSLP::serialize(out);
  }

  void load(std::istream& in) {
    DiffBase::load(in);
    TSampledSLP::load(in);
  }

 private:
  void buildSampledSLP(const sdsl::int_vector<>& da,
                       uint32_t block_size,
                       float storing_factor,
                       grammar::Chunks<>& cslp_docs_out,
                       const grammar::SLP<>* cached_cnf);
};

// RePair the DA into a CNF SLP for the sampled tree. bs/sf/variant-independent,
// so it can be cached once per collection and reused across all DGCDA cells.
inline grammar::SLP<> BuildDiffCnfSlp(const sdsl::int_vector<>& da) {
  const auto n = da.size();
  std::vector<int> da_vec(n);
  for (std::size_t i = 0; i < n; ++i)
    da_vec[i] = static_cast<int>(da[i]);
  grammar::SLP<> slp_cnf;
  grammar::RePairEncoder<true> encoder;
  auto wrapper = grammar::BuildSLPWrapper(slp_cnf);
  encoder.Encode(da_vec.begin(), da_vec.end(), wrapper);
  return slp_cnf;
}

//~~~~~~~


template <typename TSLP,
          typename TSampledSLP,
          typename TRoots,
          typename TSpanSums,
          typename TSamples,
          typename TSampleRootsPos,
          typename TBV>
void DifferentialLightSLP<TSLP, TSampledSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>::buildSampledSLP(
    const sdsl::int_vector<>& da,
    uint32_t block_size,
    float storing_factor,
    grammar::Chunks<>& cslp_docs_out,
    const grammar::SLP<>* cached_cnf) {
  // Use the cached CNF SLP if provided; otherwise RePair the DA here.
  grammar::SLP<> built_cnf;
  if (!cached_cnf) {
    built_cnf = BuildDiffCnfSlp(da);
    cached_cnf = &built_cnf;
  }

  grammar::CombinedSLP<> cslp(*cached_cnf);
  grammar::AddSet add_set(cslp_docs_out);
  cslp.Compute(block_size,
               add_set,
               add_set,
               grammar::MustBeSampled<grammar::Chunks<>>(grammar::AreChildrenTooBig(cslp_docs_out, storing_factor)));

  // Assign TSampledSLP base; rank/select pointers temporarily dangling — fixed on load from cache
  static_cast<TSampledSLP&>(*this) = static_cast<const TSampledSLP&>(cslp);
}

//~~~~~~~


// Per-field size breakdown for benchmark reporting. Delegates the DiffBase part to
// DifferentialSLP's collectSizes and adds the sampled-SLP component.
template <typename TSLP,
          typename TSampledSLP,
          typename TRoots,
          typename TSpanSums,
          typename TSamples,
          typename TSampleRootsPos,
          typename TBV>
void collectSizes(SizeReport& out,
                  const DifferentialLightSLP<TSLP, TSampledSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& slp,
                  const std::string& prefix = "") {
  using DiffBase = DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>;
  collectSizes(out, static_cast<const DiffBase&>(slp), prefix);
  append(out, prefix + "sampled_slp", sdsl::size_in_bytes(static_cast<const TSampledSLP&>(slp)));
}

//~~~~~~~


// Dedicated overload — template deduction does not cross derived→base boundaries,
// so explicitly upcast and forward to the DifferentialSLP overload.
template <typename TSLP,
          typename TSampledSLP,
          typename TRoots,
          typename TSpanSums,
          typename TSamples,
          typename TSampleRootsPos,
          typename TBV,
          typename Report>
void ExpandSLP(const DifferentialLightSLP<TSLP, TSampledSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& slp,
               std::size_t bp,
               std::size_t ep,
               Report& report) {
  ExpandSLP(static_cast<const DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>&>(slp),
            bp,
            ep,
            report);
}

//~~~~~~~


template <typename TSLP,
          typename TSampledSLP,
          typename TRoots,
          typename TSpanSums,
          typename TSamples,
          typename TSampleRootsPos,
          typename TBV>
void construct(DifferentialLightSLP<TSLP, TSampledSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& t_dslp,
               Config& t_config,
               uint32_t block_size,
               float storing_factor) {
  using namespace conf;
  using DLSLP = DifferentialLightSLP<TSLP, TSampledSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>;

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  // Both RePair passes are bs/sf-independent, so cache them once per collection and
  // reuse across all DGCDA cells:
  //  (1) the diff base grammar (shared by container variants with the same TSLP), and
  //  (2) the sampled-tree CNF SLP (a plain grammar::SLP<>, shared by all variants).
  const std::string slp_name = t_config.keys[kDGCDA][kSLP].get<std::string>();
  const std::string key_cnf = slp_name + "_cnf";
  grammar::SLP<> cnf;
  if (sdsl::cache_file_exists<grammar::SLP<>>(key_cnf, t_config)) {
    sdsl::load_from_cache(cnf, key_cnf, t_config, true);
  } else {
    cnf = BuildDiffCnfSlp(da);
    sdsl::store_to_cache(cnf, key_cnf, t_config, true);
  }

  grammar::Chunks<> cslp_docs;
  // Load-or-build the diff base grammar (skips RePair on cache hit), finish the
  // bs-dependent sampling, build the sampled tree from the cached CNF, then convert
  // into t_dslp to bit-compress the base.
  DLSLP tmp;
  auto compact_seq = LoadOrBuildDiffGrammar<TSLP>(tmp, da, t_config, slp_name + "_grammar");
  tmp.FinishCompute(block_size, compact_seq);
  tmp.BuildSampled(da, block_size, storing_factor, cslp_docs, &cnf);
  auto bit_compress = [](auto& v) {
    if constexpr (std::is_same_v<std::decay_t<decltype(v)>, sdsl::int_vector<>>)
      sdsl::util::bit_compress(v);
  };
  t_dslp = DLSLP(tmp, bit_compress, bit_compress);

  const std::string key_prefix = std::format("{}-{}_", block_size, storing_factor);
  const std::string key_docs = key_prefix + t_config.keys[kDGCDA][kDocs].get<std::string>();

  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);

  const std::string key_slp = key_prefix + t_config.keys[kDGCDA][kSLP].get<std::string>();
  sdsl::store_to_cache(t_dslp, key_slp, t_config, true);
}

}  // namespace dret
