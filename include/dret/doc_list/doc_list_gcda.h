//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/21/26.
//

#pragma once

#include <stdexcept>
#include <type_traits>

#include <sdsl/construct_sa.hpp>

#include <grammar/re_pair.h>
#include <grammar/sampled_slp.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "sr-index/r_index.h"

#include "dret/slp/combined_slp_with_unit_cover.h"
#include "dret/slp/compact_bp_slp.h"
#include "dret/slp/compact_louds_slp.h"
#include "dret/slp/differential_light_slp.h"
#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_sampled_tree_base.h"
#include "dret/index_base.h"
#include "dret/slp/slp_components.h"
#include "dret/slp/slp_tools.h"

namespace dret {

namespace gcda {

class MergeSetsBinaryTreeFunctor;

// SLP-variant trait GCDAVariantTraits<TSLP>: forward-declared here so
// DocListIdxGCDA::loadInner() can reach ::kKey via dependent-name lookup.
// The full definition (including constructSLP) lives below the
// construct(SLPType&, ...) forward decls so its body sees them via ordinary
// lookup at template-definition time.
template <typename TSLP>
struct GCDAVariantTraits;

// Key of the document lists of a cell's sampled nodes in a codec: bit-packed
// (plain) or grammar-compressed (Re-Pair). The raw lists every codec is built
// from are cached under the plain key in grammar::Chunks<>'s encoding.
template <typename TSets>
std::string NodeDocListsKey(const JSON& t_keys, const SampledTreeCell& t_cell) {
  constexpr bool kPlain = std::is_same_v<TSets, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>> ||
                          std::is_same_v<TSets, grammar::Chunks<>>;
  return CellKey(t_keys, kPlain ? conf::kDaNodeDocListsPlain : conf::kDaNodeDocListsRP, t_cell);
}

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::LightSLP<grammar::BasicSLP<sdsl::int_vector<>>,
                                            grammar::SampledSLP<>,
                                            grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TSLPSets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                                true,
                                                grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TMergeSets = MergeSetsBinaryTreeFunctor>
class DocListIdxGCDA : public DLSampledTreeScheme<TMergeSets, typename TAlphabet::string_type>,
                       public IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width> {
 public:
  using SchemeBase = DLSampledTreeScheme<TMergeSets, typename TAlphabet::string_type>;
  using StorageBase = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using typename SchemeBase::size_type;

  explicit DocListIdxGCDA(const TStorage& t_storage, uint32_t t_block_size = 512, float t_storing_factor = 4)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_storage),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor) {}

  DocListIdxGCDA(const TStorage& t_storage,
                 const TCountIdx& t_count_idx,
                 uint32_t t_block_size = 512,
                 float t_storing_factor = 4)
      : SchemeBase(TMergeSets()),
        StorageBase(t_storage),
        count_idx_(t_count_idx),
        block_size_(t_block_size),
        storing_factor_(t_storing_factor) {}

  DocListIdxGCDA() = default;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return count_idx_.serialize(out, child, "count_idx") + slp_->serialize(out, child, "slp")
           + slp_sets_->serialize(out, child, "slp_docs");
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    if (slp_) collectSizes(r, *slp_, "slp_");
    if (slp_sets_) append(r, "slp_sets", sdsl::size_in_bytes(*slp_sets_));
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    return r;
  }

  const TCountIdx& count_idx() const {
    return count_idx_;
  }

  const uint32_t& block_size() const {
    return block_size_;
  }

  const float& storing_factor() const {
    return storing_factor_;
  }

 protected:
  std::pair<std::size_t, std::size_t> count(const typename SchemeBase::TPattern& t_pattern) const override {
    return count_idx_.Count(t_pattern);
  }

  std::pair<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> computeCover(
      std::size_t t_sp,
      std::size_t t_ep) const override {
    std::vector<std::size_t> nodes;
    auto report = [&nodes](const auto& _value) {
      nodes.emplace_back(_value);
    };
    auto range = grammar::ComputeCoverFromBottom(*slp_, t_sp, t_ep, report);
    return {std::move(range), std::move(nodes)};
  }

  void getDocs(std::size_t t_sp,
               std::size_t t_ep,
               const std::function<void(typename SchemeBase::TDocId)>& t_report) const override {
    ExpandSLP(*slp_, t_sp, t_ep, t_report);
  }

  std::vector<typename SchemeBase::TDocId> getDocSet(std::size_t t_i) const override {
    auto v = (*slp_sets_)[t_i];
    return {v.begin(), v.end()};
  }

  void loadInner(typename StorageBase::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit(
        [this](auto&& tt_source) {
          count_idx_.load(tt_source.get());
        },
        t_source);

    const SampledTreeCell cell{block_size_, storing_factor_};
    slp_ = this->template loadComponentsPtr<TSLP>(CacheComponents(TSLP{}, t_keys, cell), t_source);
    slp_sets_ = this->template loadItemPtr<TSLPSets>(NodeDocListsKey<TSLPSets>(t_keys, cell), t_source, true);
  }

  TCountIdx count_idx_;
  const TSLP* slp_ = nullptr;
  const TSLPSets* slp_sets_ = nullptr;

  uint32_t block_size_ = 512;
  float storing_factor_ = 4;
};

//~~~~~~~


template <bool kExpand, typename TChunks>
void construct(grammar::GCChunks<grammar::SLP<>, kExpand, TChunks>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename TSLP, bool kExpand>
void construct(grammar::GCChunks<TSLP, kExpand, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor);

// Plain (bit-packed, not grammar-compressed) document sets: the sorted lists
// themselves, one chunk per sampled node, as PDL's plain codec stores them.
inline void construct(grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>& t_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor);

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

// Compact-grammar TSLP variants (added in Phase A.0.3 upstream; the dret-side
// shims at "dret/slp/compact_bp_slp.h" / "dret/slp/compact_louds_slp.h" re-export them as
// dret::CompactBPSLP / dret::CompactLOUDSSLP). These forward declarations are
// required so the unqualified `construct(slp, ...)` call inside the
// `construct(DocListIdxGCDA&, ...)` template below resolves at the
// template-definition point — ordinary lookup is fixed there, and ADL on
// `grammar::CompactBPSLP` searches `namespace grammar` only, never
// `dret::gcda` where the actual overload definitions live.
template <typename... Ts>
void construct(grammar::CompactBPSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

template <typename... Ts>
void construct(grammar::CompactLOUDSSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

// Phase B: CombinedSLP-with-unit-cover variant. Same signature pattern as
// the compact-grammar variants above; same forward-decl placement reason
// (visibility for the dispatch template's unqualified `construct(slp, ...)`
// call below).
template <typename... Ts>
void construct(grammar::CombinedSLPWithUnitCover<Ts...>& t_wrapped,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor);

//~~~~~~~


// SLP-variant trait full definition. Lives below the construct(SLPType&, ...)
// forward declarations so its constructSLP body sees them via ordinary lookup
// at template-definition time. The primary template covers every
// GCDA-family SLP (grammar::LightSLP, CompactBPSLP, CompactLOUDSSLP,
// CombinedSLPWithUnitCover, ...): the 5-arg SLP-build path with a datafile
// parameter. The DifferentialLightSLP specialisation routes to the 4-arg
// DifferentialLightSLP::construct(). Both cache through CacheComponents.
template <typename TSLP>
struct GCDAVariantTraits {
  static void constructSLP(TSLP& slp,
                           Config& t_config,
                           const std::string& t_datafile,
                           uint32_t t_block_size,
                           float t_storing_factor) {
    construct(slp, t_config, t_datafile, t_block_size, t_storing_factor);
  }
};

template <typename... Ts>
struct GCDAVariantTraits<DifferentialLightSLP<Ts...>> {
  static void constructSLP(DifferentialLightSLP<Ts...>& slp,
                           Config& t_config,
                           const std::string& /*t_datafile_ignored*/,
                           uint32_t t_block_size,
                           float t_storing_factor) {
    // DifferentialLightSLP::construct (defined in
    // dret/slp/differential_light_slp.h) loads the DA from the config cache
    // itself, so no datafile parameter.
    construct(slp, t_config, t_block_size, t_storing_factor);
  }
};


template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TSLP,
          typename TSLPSets,
          typename TMergeSets>
void construct(DocListIdxGCDA<TStorage, TAlphabet, TCountIdx, TSLP, TSLPSets, TMergeSets>& t_index, Config& t_config) {
  using namespace conf;
  using Traits = GCDAVariantTraits<TSLP>;

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

  // ConstructDocArray writes DA only with type-hash; the no-hash check would never
  // see it, causing DA to be rebuilt on every construct() call. Match the storage.
  if (const auto key = t_config.keys[kDA].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray(t_config);
  }

  const SampledTreeCell cell{t_index.block_size(), t_index.storing_factor()};

  if (const auto components = CacheComponents(TSLP{}, t_config.keys, cell); !ComponentsExist(components, t_config)) {
    auto event = sdsl::memory_monitor::event(components.back().key);
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
    TSLP slp;
    Traits::constructSLP(slp, t_config, filepath_da, t_index.block_size(), t_index.storing_factor());
  }

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  // The lists are the same for every representation of the cell (built while
  // sampling the tree, which every representation shares).
  if (const auto key_docs = NodeDocListsKey<TSLPSets>(t_config.keys, cell);
      !sdsl::cache_file_exists<TSLPSets>(key_docs, t_config)) {
    auto event = sdsl::memory_monitor::event(key_docs);
    TSLPSets slp_sets;
    construct(slp_sets, t_config, t_index.block_size(), t_index.storing_factor());
  }

  t_index.load(t_config);
}

//~~~~~~~


// The sampled tree of a (block size, storing factor) cell over the CNF grammar
// of the DA, with the document lists of its sampled nodes. Every representation
// derives from it; it is cached as components (the CNF grammar, the tree, its
// leaves), and the raw lists under the plain lists key, in grammar::Chunks<>'s
// encoding and bit-packed.
template <typename TSLP, typename TSampledSLP, typename TLeavesContainer>
void construct(grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>& t_cslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto components = CacheComponents(t_cslp, t_config.keys, cell);
  auto event = sdsl::memory_monitor::event(ComponentFile(components.back(), t_config));

  t_cslp = grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>(LoadOrBuildDaCnfGrammar(t_config, t_datafile));
  grammar::Chunks<> cslp_docs;

  grammar::AddSet add_set(cslp_docs);
  t_cslp.Compute(t_block_size,
                 add_set,
                 add_set,
                 grammar::MustBeSampled<decltype(cslp_docs)>(grammar::AreChildrenTooBig(cslp_docs, t_storing_factor)));

  StoreComponents(t_cslp, CacheComponents(t_cslp, t_config.keys, cell), t_config);

  const auto key_docs = CellKey(t_config.keys, kDaNodeDocListsPlain, cell);
  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);
}

// The raw sampled tree of a cell (see above): loaded from its components, or built.
inline grammar::CombinedSLP<> LoadOrBuildSampledTree(Config& t_config,
                                                     const std::string& t_datafile,
                                                     uint32_t t_block_size,
                                                     float t_storing_factor) {
  grammar::CombinedSLP<> cslp;
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  if (!LoadComponents(cslp, CacheComponents(cslp, t_config.keys, cell), t_config))
    construct(cslp, t_config, t_datafile, t_block_size, t_storing_factor);
  return cslp;
}

//~~~~~~~


// The SLP over the document array, plus the compact sequence, as parsed from
// irepair's .R/.C. Cached once per collection under bs/sf-independent keys.
//
// irepair itself already ran once per collection -- every call site guards it
// on std::filesystem::exists(t_datafile + ".R"), a path with no (bs,sf) in it.
// What used to repeat was this *parse* of its output: it lived inside the
// per-cell guard in construct() below, so all 20 cells of a grid sweep redid
// it. Measured across the 2026-08 campaign that cost ~71 h, ~8.4% of stage A,
// and it is why revision-mid spent 12.85 h re-deriving one grammar 17 times.
//
// compact_seq round-trips through a bit-compressed int_vector; the values are
// grammar variable ids, so they are non-negative and bounded by the rule count.
inline void LoadOrBuildDaSlp(Config& t_config,
                             const std::string& t_datafile,
                             grammar::SLP<>& t_slp,
                             std::vector<std::size_t>& t_compact_seq) {
  using namespace conf;
  const auto key_slp = KeyName(t_config.keys, kDaGrammar);
  const auto key_seq = KeyName(t_config.keys, kDaGrammarSeq);

  if (sdsl::cache_file_exists<grammar::SLP<>>(key_slp, t_config)
      && sdsl::cache_file_exists<sdsl::int_vector<>>(key_seq, t_config)) {
    sdsl::load_from_cache(t_slp, key_slp, t_config, true);
    sdsl::int_vector<> seq;
    sdsl::load_from_cache(seq, key_seq, t_config, true);
    t_compact_seq.assign(seq.begin(), seq.end());
    return;
  }

  auto event = sdsl::memory_monitor::event(key_slp);
  {
    grammar::RePairReader<false> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(t_slp);
    auto report_compact_seq = [&t_compact_seq](const auto& _var) {
      t_compact_seq.emplace_back(_var);
    };
    re_pair_reader.Read(t_datafile, slp_wrapper, report_compact_seq);
  }

  sdsl::store_to_cache(t_slp, key_slp, t_config, true);
  sdsl::int_vector<> seq(t_compact_seq.size());
  std::copy(t_compact_seq.begin(), t_compact_seq.end(), seq.begin());
  sdsl::util::bit_compress(seq);
  sdsl::store_to_cache(seq, key_seq, t_config, true);
}

//~~~~~~~


// Sampled-cached: the DA grammar (not in CNF) with, for each leaf of the cell's
// sampled tree, its cover in that grammar.
template <typename TSLP, typename TSampledSLP, typename TChunks>
void construct(grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto cslp = LoadOrBuildSampledTree(t_config, t_datafile, t_block_size, t_storing_factor);

  // Collection-level, not cell-level: parsed once and reused by every
  // (bs,sf) cell. See LoadOrBuildDaSlp above.
  grammar::SLP<> slp;
  std::vector<std::size_t> compact_seq;
  LoadOrBuildDaSlp(t_config, t_datafile, slp, compact_seq);

  auto event = sdsl::memory_monitor::event(ComponentFile(CacheComponents(t_lslp, t_config.keys, cell).back(), t_config));
  grammar::LightSLP<> lslp;
  lslp.Compute(slp, compact_seq, cslp);

  auto bit_compress = [](auto& _v) {
    sdsl::util::bit_compress(_v);
  };
  t_lslp = grammar::LightSLP<TSLP, TSampledSLP, TChunks>(lslp, bit_compress, bit_compress, bit_compress, bit_compress);
  StoreComponents(t_lslp, CacheComponents(t_lslp, t_config.keys, cell), t_config);
}

//~~~~~~~


// Shared helper for the compact-grammar construct overloads. Both BP and
// LOUDS variants derive the compact class from the cell's raw sampled tree via
// Compute(): a compact encoding of the CNF grammar (one per collection), the
// sampled leaves numbered in it, and the shared sampled tree.
template <typename TCompactSLP>
void constructCompactCommon(TCompactSLP& t_compact,
                            Config& t_config,
                            const std::string& t_datafile,
                            uint32_t t_block_size,
                            float t_storing_factor) {
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto cslp = LoadOrBuildSampledTree(t_config, t_datafile, t_block_size, t_storing_factor);

  auto event = sdsl::memory_monitor::event(ComponentFile(CacheComponents(t_compact, t_config.keys, cell)[1], t_config));
  t_compact.Compute(cslp);
  StoreComponents(t_compact, CacheComponents(t_compact, t_config.keys, cell), t_config);
}

template <typename... Ts>
void construct(grammar::CompactBPSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  constructCompactCommon(t_compact, t_config, t_datafile, t_block_size, t_storing_factor);
}

template <typename... Ts>
void construct(grammar::CompactLOUDSSLP<Ts...>& t_compact,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  constructCompactCommon(t_compact, t_config, t_datafile, t_block_size, t_storing_factor);
}

// Sampled-ondemand: the cell's raw sampled tree converted into the wrapper's
// (possibly bit-compressed) CombinedSLP base -- applying `bit_compress` to the
// SLP rules, lengths, and leaves.
template <typename... Ts>
void construct(grammar::CombinedSLPWithUnitCover<Ts...>& t_wrapped,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using BaseCSLP = typename grammar::CombinedSLPWithUnitCover<Ts...>::Base;
  const SampledTreeCell cell{t_block_size, t_storing_factor};
  const auto raw = LoadOrBuildSampledTree(t_config, t_datafile, t_block_size, t_storing_factor);

  auto event = sdsl::memory_monitor::event(ComponentFile(CacheComponents(t_wrapped, t_config.keys, cell).back(), t_config));

  // The production SLP_Combined uses sdsl::int_vector<> containers; the action is a
  // no-op for any non-int_vector container (e.g. the default std::vector base).
  auto bit_compress = [](auto& v) {
    if constexpr (std::is_same_v<std::decay_t<decltype(v)>, sdsl::int_vector<>>)
      sdsl::util::bit_compress(v);
  };
  BaseCSLP cslp(raw, bit_compress, bit_compress, bit_compress);

  static_cast<BaseCSLP&>(t_wrapped) = cslp;
  StoreComponents(t_wrapped, CacheComponents(t_wrapped, t_config.keys, cell), t_config);
}

//~~~~~~~


// The raw lists of a cell, cached while its sampled tree was built.
inline grammar::Chunks<> LoadRawNodeDocLists(Config& t_config, uint32_t t_block_size, float t_storing_factor) {
  const auto key_docs = NodeDocListsKey<grammar::Chunks<>>(t_config.keys, {t_block_size, t_storing_factor});
  grammar::Chunks<> sets;
  if (!sdsl::load_from_cache(sets, key_docs, t_config, true))
    throw std::runtime_error("GCDA document lists: no cached lists under " + key_docs +
                             "; build the sampled tree for this (block size, storing factor) first");
  return sets;
}

template <bool kExpand, typename TChunks>
void construct(grammar::GCChunks<grammar::SLP<>, kExpand, TChunks>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor) {
  const auto slp_sets = LoadRawNodeDocLists(t_config, t_block_size, t_storing_factor);
  const auto& objs = slp_sets.GetObjects();
  grammar::RePairEncoder<false> encoder_nslp;
  t_slp_sets.Compute(objs.begin(), objs.end(), slp_sets, encoder_nslp);
  sdsl::store_to_cache(t_slp_sets, NodeDocListsKey<std::decay_t<decltype(t_slp_sets)>>(t_config.keys,
                                                                                        {t_block_size, t_storing_factor}),
                       t_config, true);
}

//~~~~~~~


template <typename TSLP, bool kExpand>
void construct(grammar::GCChunks<TSLP, kExpand, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor) {
  const auto key_docs = NodeDocListsKey<std::decay_t<decltype(t_slp_sets)>>(t_config.keys,
                                                                            {t_block_size, t_storing_factor});

  grammar::GCChunks<grammar::SLP<>> slp_sets;
  if (!sdsl::cache_file_exists<decltype(slp_sets)>(key_docs, t_config)) {
    construct(slp_sets, t_config, t_block_size, t_storing_factor);
  } else {
    sdsl::load_from_cache(slp_sets, key_docs, t_config, true);
  }

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };

  t_slp_sets =
      std::remove_reference_t<decltype(t_slp_sets)>(slp_sets, bit_compress, bit_compress, bit_compress, bit_compress);
  sdsl::store_to_cache(t_slp_sets, key_docs, t_config, true);
}

//~~~~~~~


// Plain (bit-packed, not grammar-compressed) document lists: the sorted lists
// themselves, one chunk per sampled node, as PDL's plain codec stores them.
// Cached as a side effect of building the sampled tree; GCDA's construct
// reaches this only when the bit-packed copy is missing.
inline void construct(grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>& t_sets,
                      Config& t_config,
                      uint32_t t_block_size,
                      float t_storing_factor) {
  const auto sets = LoadRawNodeDocLists(t_config, t_block_size, t_storing_factor);
  auto bit_compress = [](sdsl::int_vector<>& v) { sdsl::util::bit_compress(v); };
  t_sets = grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>(sets, bit_compress, bit_compress);
  sdsl::store_to_cache(t_sets, NodeDocListsKey<grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>(
                                   t_config.keys, {t_block_size, t_storing_factor}),
                       t_config, true);
}

//~~~~~~~


class MergeSetsBinaryTreeFunctor {
 public:
  MergeSetsBinaryTreeFunctor() = default;

  template <typename TII, typename TSets, typename TResult>
  void operator()(TII _first, TII _last, const TSets& _sets, TResult& _result) const {
    auto default_set_union = [](auto _first1, auto _last1, auto _first2, auto _last2, auto _result) -> auto {
      return std::set_union(_first1, _last1, _first2, _last2, _result);
    };

    (*this)(_first, _last, _sets, _result, default_set_union);
  }

  template <typename _II, typename _Sets, typename _Result, typename _SetUnion>
  void operator()(_II _first, _II _last, const _Sets& _sets, _Result& _result, const _SetUnion& _set_union) const {
    _Result tmp_merge;

    auto merge_tmp = [&tmp_merge, &_set_union](const auto& set1, const auto& set2) {
      tmp_merge.resize(set1.size() + set2.size());
      auto last_it = _set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), tmp_merge.begin());
      tmp_merge.resize(last_it - tmp_merge.begin());
    };

    auto length = std::distance(_first, _last);
    if (length == 1) {
      merge_tmp(_result, _sets(*_first));

      _result.swap(tmp_merge);

      return;
    }

    std::vector<std::pair<uint8_t, _Result>> part_results = {{1, {}}};
    part_results.front().second.swap(_result);

    while (part_results.size() != 1 || _first != _last) {
      std::size_t size;
      while ((size = part_results.size()) > 1
             && (part_results[size - 1].first == part_results[size - 2].first || _first == _last)) {
        merge_tmp(part_results[size - 1].second, part_results[size - 2].second);

        part_results[size - 2].second.swap(tmp_merge);
        ++part_results[size - 2].first;
        part_results.pop_back();
      }

      if (_first != _last) {
        auto next = _first + 1;
        if (next == _last) {
          merge_tmp(part_results.back().second, _sets(*_first));

          part_results.back().second.swap(tmp_merge);
          ++_first;
        } else {
          merge_tmp(_sets(*_first), _sets(*next));

          part_results.emplace_back(1, std::move(tmp_merge));
          _first += 2;
        }
      }
    }

    _result.swap(part_results.front().second);
  }
};

}  // namespace gcda

//~~~~~~~

// Backward-compatible alias. dret::dgcda::DocListIdxDGCDA locks TSLP to
// DifferentialLightSLP<> and forwards every other template parameter and
// default to gcda::DocListIdxGCDA. Existing consumers that wrote
// dret::dgcda::DocListIdxDGCDA<...> keep working unchanged; ADL on a
// DocListIdxGCDA<..., DifferentialLightSLP<>, ...> argument resolves to the
// unified gcda::construct(DocListIdxGCDA&, Config&) above.
namespace dgcda {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = DifferentialLightSLP<>,
          typename TSLPSets = grammar::GCChunks<grammar::BasicSLP<sdsl::int_vector<>>,
                                                  true,
                                                  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>,
          typename TMergeSets = gcda::MergeSetsBinaryTreeFunctor>
using DocListIdxDGCDA = gcda::DocListIdxGCDA<TStorage, TAlphabet, TCountIdx,
                                              TSLP, TSLPSets, TMergeSets>;

}  // namespace dgcda

}  // namespace dret
