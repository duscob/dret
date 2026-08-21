//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 3/21/26.
//

#pragma once

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

    auto key_prefix = std::format("{}-{}_", block_size_, storing_factor_);
    constexpr auto kVariantKey = GCDAVariantTraits<TSLP>::kKey;
    slp_ = this->template loadItemPtr<TSLP>(
        key_prefix + t_keys[kVariantKey][kSLP].template get<std::string>(), t_source, true);
    slp_sets_ = this->template loadItemPtr<TSLPSets>(
        key_prefix + t_keys[kVariantKey][kDocs].template get<std::string>(), t_source, true);
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
// CombinedSLPWithUnitCover, ...): cache keys under conf::kGCDA and the
// 5-arg SLP-build path with a datafile parameter. The DifferentialLightSLP
// specialisation routes to the DGCDA cache namespace and the 4-arg
// DifferentialLightSLP::construct().
template <typename TSLP>
struct GCDAVariantTraits {
  static constexpr std::string_view kKey = conf::kGCDA;

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
  static constexpr std::string_view kKey = conf::kDGCDA;

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

  const auto key_prefix = std::format("{}-{}_", t_index.block_size(), t_index.storing_factor());

  if (const auto key = key_prefix + t_config.keys[Traits::kKey][kSLP].template get<std::string>();
      !sdsl::cache_file_exists<TSLP>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
    TSLP slp;
    Traits::constructSLP(slp, t_config, filepath_da, t_index.block_size(), t_index.storing_factor());
  }

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  if (const auto key_docs = key_prefix + t_config.keys[Traits::kKey][kDocs].template get<std::string>();
      !sdsl::cache_file_exists<TSLPSets>(key_docs, t_config)) {
    auto event = sdsl::memory_monitor::event(key_docs);
    if constexpr (Traits::kKey == conf::kDGCDA) {
      // DifferentialLightSLP::construct stored a plain grammar::Chunks<> at
      // key_docs as a side effect; build the GCChunks-of-Chunks from it
      // directly. The GCDA dispatched path (else branch) hardcodes kGCDA in
      // its cache lookups and so cannot serve the DGCDA key prefix.
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
    } else {
      TSLPSets slp_sets;
      construct(slp_sets, t_config, t_index.block_size(), t_index.storing_factor());
    }
  }

  t_index.load(t_config);
}

//~~~~~~~


template <typename TSLP, typename TSampledSLP, typename TLeavesContainer>
void construct(grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>& t_cslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;

  // Grammar compress data file using RePair
  if (!std::filesystem::exists(t_datafile + ".R") && REPAIR_EXE) {
    const auto filename = std::filesystem::path(t_datafile).filename().string();
    auto event = sdsl::memory_monitor::event("RePair-" + filename);
    RunRePair(t_datafile);
    t_config.file_map[filename + ".R"] = t_datafile + ".R";
    t_config.file_map[filename + ".C"] = t_datafile + ".C";
  }
  CheckRePairGrammar(t_datafile);

  const auto key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);

  auto key_slp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>>(key_slp, t_config));

  grammar::SLP<> slp;
  {
    grammar::RePairReader<true> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(slp);
    re_pair_reader.Read(t_datafile, slp_wrapper);
  }

  t_cslp = grammar::CombinedSLP<TSLP, TSampledSLP, TLeavesContainer>(slp);
  grammar::Chunks<> cslp_docs;

  grammar::AddSet add_set(cslp_docs);
  t_cslp.Compute(t_block_size,
                 add_set,
                 add_set,
                 grammar::MustBeSampled<decltype(cslp_docs)>(grammar::AreChildrenTooBig(cslp_docs, t_storing_factor)));

  sdsl::store_to_cache(t_cslp, key_slp, t_config, true);

  auto key_docs = key_prefix + t_config.keys[kGCDA][kDocs].get<std::string>();
  sdsl::store_to_cache(cslp_docs, key_docs, t_config, true);

  auto bit_compress = [](sdsl::int_vector<>& _v) {
    sdsl::util::bit_compress(_v);
  };
  grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>> cslp_docs_c(cslp_docs, bit_compress, bit_compress);
  sdsl::store_to_cache(cslp_docs_c, key_docs, t_config, true);
}

//~~~~~~~


// The SLP over the document array, plus the compact sequence, as parsed from
// irepair's .R/.C. Cached once per collection under bs/sf-independent keys.
//
// irepair itself already ran once per collection -- every call site guards it
// on std::filesystem::exists(t_datafile + ".R"), a path with no (bs,sf) in it.
// What used to repeat was this *parse* of its output: it lived inside the
// "<bs>-<sf>_gcda_slp" guard in construct() below, so all 20 cells of a grid
// sweep redid it. Measured across the 2026-08 campaign that cost ~71 h, ~8.4%
// of stage A, and it is why revision-mid spent 12.85 h re-deriving one grammar
// 17 times. DGCDA never had the problem -- see dgcda_slp_grammar.
//
// compact_seq round-trips through a bit-compressed int_vector; the values are
// grammar variable ids, so they are non-negative and bounded by the rule count.
inline void LoadOrBuildDaSlp(Config& t_config,
                             const std::string& t_datafile,
                             grammar::SLP<>& t_slp,
                             std::vector<std::size_t>& t_compact_seq) {
  using namespace conf;
  const auto key_slp = t_config.keys[kGCDA][kSLPGrammar].get<std::string>();
  const auto key_seq = t_config.keys[kGCDA][kSLPCompactSeq].get<std::string>();

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


template <typename TSLP, typename TSampledSLP, typename TChunks>
void construct(grammar::LightSLP<TSLP, TSampledSLP, TChunks>& t_lslp,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;

  std::string key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);

  auto key_lslp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
  grammar::LightSLP<> lslp;

  if (!sdsl::cache_file_exists<decltype(lslp)>(key_lslp, t_config)) {
    // Construct Light SLP on datafile
    auto event = sdsl::memory_monitor::event(sdsl::cache_file_name<grammar::LightSLP<>>(key_lslp, t_config));

    grammar::CombinedSLP<> cslp;
    if (const auto key = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();
        !sdsl::cache_file_exists<decltype(cslp)>(key, t_config)) {
      auto event_cslp = sdsl::memory_monitor::event(sdsl::cache_file_name<decltype(cslp)>(key, t_config));
      construct(cslp, t_config, t_datafile, t_block_size, t_storing_factor);
    } else {
      sdsl::load_from_cache(cslp, key, t_config, true);
    }

    // Collection-level, not cell-level: parsed once and reused by every
    // (bs,sf) cell. See LoadOrBuildDaSlp above.
    grammar::SLP<> slp;
    std::vector<std::size_t> compact_seq;
    LoadOrBuildDaSlp(t_config, t_datafile, slp, compact_seq);

    lslp.Compute(slp, compact_seq, cslp);
    sdsl::store_to_cache(lslp, key_lslp, t_config, true);
  } else {
    sdsl::load_from_cache(lslp, key_lslp, t_config, true);
  }

  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::LightSLP<TSLP, TSampledSLP, TChunks>>(key_lslp, t_config));

  // Construct Light SLP Basic on DA
  auto bit_compress = [](auto& _v) {
    sdsl::util::bit_compress(_v);
  };
  t_lslp = grammar::LightSLP<TSLP, TSampledSLP, TChunks>(lslp, bit_compress, bit_compress, bit_compress, bit_compress);
  sdsl::store_to_cache(t_lslp, key_lslp, t_config, true);
}

//~~~~~~~


// Shared helper for the compact-grammar construct overloads. Both BP and
// LOUDS variants follow the same pattern: ensure a default
// grammar::CombinedSLP<> is present in the cache (building it via the
// existing construct(CombinedSLP&,...) overload if needed — that path
// also writes the gcda_docs cache, which the subsequent
// construct(GCChunks<>,...) consumes), then derive the compact class from
// the CSLP via Compute(), and store it under the same logical SLP key
// (the type-hash in store_to_cache disambiguates the on-disk file name
// across variants).
template <typename TCompactSLP>
void constructCompactCommon(TCompactSLP& t_compact,
                            Config& t_config,
                            const std::string& t_datafile,
                            uint32_t t_block_size,
                            float t_storing_factor) {
  using namespace conf;

  std::string key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);
  auto key_slp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();

  grammar::CombinedSLP<> cslp;
  if (!sdsl::cache_file_exists<grammar::CombinedSLP<>>(key_slp, t_config)) {
    // Lexical scope is dret::gcda; this resolves to the existing
    // construct(grammar::CombinedSLP&,...) overload above, which also
    // writes the gcda_docs cache.
    construct(cslp, t_config, t_datafile, t_block_size, t_storing_factor);
  } else {
    sdsl::load_from_cache(cslp, key_slp, t_config, true);
  }

  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<TCompactSLP>(key_slp, t_config));

  t_compact.Compute(cslp);
  sdsl::store_to_cache(t_compact, key_slp, t_config, true);
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

// Phase B: CombinedSLP-with-unit-cover construct. Builds/loads a default
// std::vector-backed `grammar::CombinedSLP<>` via the existing CSLP construct
// overload (which writes the `gcda_docs` cache), converts it into the wrapper's
// (possibly bit-compressed) CombinedSLP base — applying `bit_compress` to the
// SLP rules, lengths, and leaves — slices that into the wrapper, then stores the
// wrapper under the same logical SLP key. Type-hashing keeps the wrapper's
// on-disk file distinct from the raw CSLP's.
template <typename... Ts>
void construct(grammar::CombinedSLPWithUnitCover<Ts...>& t_wrapped,
               Config& t_config,
               const std::string& t_datafile,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;
  using BaseCSLP = typename grammar::CombinedSLPWithUnitCover<Ts...>::Base;
  using RawCSLP = grammar::CombinedSLP<>;  // std::vector intermediate built by construct(CombinedSLP&,...)

  std::string key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);
  auto key_slp = key_prefix + t_config.keys[kGCDA][kSLP].get<std::string>();

  RawCSLP raw;
  if (!sdsl::cache_file_exists<RawCSLP>(key_slp, t_config)) {
    construct(raw, t_config, t_datafile, t_block_size, t_storing_factor);
  } else {
    sdsl::load_from_cache(raw, key_slp, t_config, true);
  }

  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::CombinedSLPWithUnitCover<Ts...>>(key_slp, t_config));

  // Convert raw (std::vector) -> BaseCSLP, bit-compressing rules, lengths, leaves.
  // The production SLP_Combined uses sdsl::int_vector<> containers; the action is a
  // no-op for any non-int_vector container (e.g. the default std::vector base).
  auto bit_compress = [](auto& v) {
    if constexpr (std::is_same_v<std::decay_t<decltype(v)>, sdsl::int_vector<>>)
      sdsl::util::bit_compress(v);
  };
  BaseCSLP cslp(raw, bit_compress, bit_compress, bit_compress);

  static_cast<BaseCSLP&>(t_wrapped) = cslp;
  sdsl::store_to_cache(t_wrapped, key_slp, t_config, true);
}

//~~~~~~~


template <bool kExpand, typename TChunks>
void construct(grammar::GCChunks<grammar::SLP<>, kExpand, TChunks>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;
  const auto key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);
  auto key_docs = key_prefix + t_config.keys[kGCDA][kDocs].get<std::string>();

  grammar::Chunks<> slp_sets;
  sdsl::load_from_cache(slp_sets, key_docs, t_config, true);

  const auto& objs = slp_sets.GetObjects();
  grammar::RePairEncoder<false> encoder_nslp;
  t_slp_sets.Compute(objs.begin(), objs.end(), slp_sets, encoder_nslp);
  sdsl::store_to_cache(t_slp_sets, key_docs, t_config, true);
}

//~~~~~~~


template <typename TSLP, bool kExpand>
void construct(grammar::GCChunks<TSLP, kExpand, grammar::Chunks<sdsl::int_vector<>, sdsl::int_vector<>>>& t_slp_sets,
               Config& t_config,
               uint32_t t_block_size,
               float t_storing_factor) {
  using namespace conf;
  const auto key_prefix = std::format("{}-{}_", t_block_size, t_storing_factor);
  auto key_docs = key_prefix + t_config.keys[kGCDA][kDocs].get<std::string>();

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
