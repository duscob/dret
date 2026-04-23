//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/22/26.
//

#pragma once

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "config.h"
#include "differential_slp.h"

namespace dret {

template <typename TSLP = grammar::SLP<>,
          typename TSampledSLP = grammar::SampledSLP<>,
          typename TIntContainer = sdsl::int_vector<>,
          typename TBV = sdsl::sd_vector<>>
class DifferentialLightSLP : public DifferentialSLP<TSLP, TIntContainer, TBV>, public TSampledSLP {
  using DiffBase = DifferentialSLP<TSLP, TIntContainer, TBV>;

 public:
  using size_type = std::size_t;

  void Compute(const sdsl::int_vector<>& da,
               uint32_t block_size,
               float storing_factor,
               grammar::Chunks<>& cslp_docs_out) {
    DiffBase::Compute(da, block_size);
    buildSampledSLP(da, block_size, storing_factor, cslp_docs_out);
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
                       grammar::Chunks<>& cslp_docs_out);
};

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV>
void DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>::buildSampledSLP(const sdsl::int_vector<>& da,
                                                                                  uint32_t block_size,
                                                                                  float storing_factor,
                                                                                  grammar::Chunks<>& cslp_docs_out) {
  const auto n = da.size();

  std::vector<int> da_vec(n);
  for (std::size_t i = 0; i < n; ++i)
    da_vec[i] = static_cast<int>(da[i]);

  grammar::SLP<> slp_cnf;
  {
    grammar::RePairEncoder<true> encoder;
    auto wrapper = grammar::BuildSLPWrapper(slp_cnf);
    encoder.Encode(da_vec.begin(), da_vec.end(), wrapper);
  }

  grammar::CombinedSLP<> cslp(slp_cnf);
  grammar::AddSet add_set(cslp_docs_out);
  cslp.Compute(block_size,
               add_set,
               add_set,
               grammar::MustBeSampled<grammar::Chunks<>>(grammar::AreChildrenTooBig(cslp_docs_out, storing_factor)));

  // Assign TSampledSLP base; rank/select pointers temporarily dangling — fixed on load from cache
  static_cast<TSampledSLP&>(*this) = static_cast<const TSampledSLP&>(cslp);
}

//~~~~~~~


// Dedicated overload — template deduction does not cross derived→base boundaries,
// so explicitly upcast and forward to the DifferentialSLP overload.
template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV, typename Report>
void ExpandSLP(const DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>& slp,
               std::size_t bp,
               std::size_t ep,
               Report& report) {
  ExpandSLP(static_cast<const DifferentialSLP<TSLP, TIntContainer, TBV>&>(slp), bp, ep, report);
}

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TIntContainer, typename TBV>
void construct(DifferentialLightSLP<TSLP, TSampledSLP, TIntContainer, TBV>& t_dslp,
               Config& t_config,
               uint32_t block_size,
               float storing_factor) {
  using namespace conf;

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  grammar::Chunks<> cslp_docs;
  t_dslp.Compute(da, block_size, storing_factor, cslp_docs);

  const std::string key_prefix = std::format("{}-{}_", block_size, storing_factor);
  const std::string key_docs = key_prefix + t_config.keys[kDGCDA][kDocs].get<std::string>();

  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);
  auto bit_compress = [](sdsl::int_vector<>& v) {
    sdsl::util::bit_compress(v);
  };
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);

  const std::string key_slp = key_prefix + t_config.keys[kDGCDA][kSLP].get<std::string>();
  sdsl::store_to_cache(t_dslp, key_slp, t_config, true);
}

}  // namespace dret
