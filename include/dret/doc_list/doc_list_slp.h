//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/26/26.
//
// Phase C: non-sampled GCDA index. Drives an SA-range search through a
// bare `grammar::SLP<>` (no sampled tree, no precomputed covers, no
// `GCChunks` doc-set storage), expanding the matched range and de-duping
// the resulting doc-id stream.
//
// Compared to `dret::gcda::DocListIdxGCDA`, this class:
//   - Has no `block_size` / `storing_factor` knobs (those control the
//     sampled tree, which is absent here).
//   - Uses `grammar::ComputeSpanCover` to split `[sp, ep)` into a small
//     set of cover variables (each spans an exact subrange), then
//     expands each via `grammar::ExpandSLPForward` to enumerate
//     terminals = doc-ids.
//   - Sorts and deduplicates the doc-id stream before reporting (mirrors
//     the dedup `DLSampledTreeScheme::Search` does after `getDocs`).
//
// `Search` cost is O(|cover| * log n) for the cover plus O((ep - sp) *
// height) for the expansions plus O(k log k) for sort/dedup with `k =
// ep - sp`. There is no `getDocSet` shortcut because there are no
// precomputed sets — the whole [sp, ep) is walked.
//

#pragma once

#include <algorithm>
#include <filesystem>
#include <format>
#include <functional>
#include <iterator>
#include <vector>

#include <sdsl/construct_sa.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/memory_management.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "dret/construct_base.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/index_base.h"
#include "dret/slp/differential_slp.h"

namespace dret {

// Forward decl of the plain bare-SLP RePair build (defined at end of file).
template <typename TVarsContainer, typename TLengthsContainer>
void construct(grammar::SLP<TVarsContainer, TLengthsContainer>& t_slp,
               Config& t_config, const std::string& t_datafile);

// Fixed differential-reconstruction sample granularity for the non-sampled
// differential SLP (bare-diff). Differential decode needs *some* anchors, so
// this is an internal constant — not an exposed bs/sf axis.
inline constexpr std::uint32_t kDiffBlockSize = 512;

//~~~~~~~  SLP-NS range expansion, dispatched by SLP type  ~~~~~~~

// Plain grammar::SLP: split [sp, ep) into a span cover, expand each variable.
template <typename TVarsContainer, typename TLengthsContainer, typename TReport>
void SlpNsExpandRange(const grammar::SLP<TVarsContainer, TLengthsContainer>& slp,
                      std::size_t sp, std::size_t ep, TReport add) {
  std::vector<typename grammar::SLP<TVarsContainer, TLengthsContainer>::VariableType> cover;
  grammar::ComputeSpanCover(slp, sp, ep, std::back_inserter(cover));
  for (auto var : cover) {
    auto length = slp.SpanLength(var);
    grammar::ExpandSLPForward(slp.GetRules(), slp.Sigma(), var, length, add);
  }
}

// Base differential SLP: range-decompress [sp, ep) directly (no sampled tree).
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples,
          typename TSampleRootsPos, typename TBV, typename TReport>
void SlpNsExpandRange(const DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>& slp,
                      std::size_t sp, std::size_t ep, TReport add) {
  ExpandSLP(slp, sp, ep, add);
}

//~~~~~~~  SLP-NS cache key + build, dispatched by SLP type  ~~~~~~~

template <typename TSLP>
struct SlpNsTraits;  // primary left undefined — only the two SLP kinds below

// Plain bare grammar::SLP — RePair build under kSLPNS (unchanged behaviour).
template <typename TVarsContainer, typename TLengthsContainer>
struct SlpNsTraits<grammar::SLP<TVarsContainer, TLengthsContainer>> {
  using TSLP = grammar::SLP<TVarsContainer, TLengthsContainer>;
  static constexpr std::string_view kKey = conf::kSLPNS;
  static void build(Config& t_config) {
    const auto key = t_config.keys[conf::kSLPNS].get<std::string>();
    if (sdsl::cache_file_exists<TSLP>(key, t_config)) return;
    auto event = sdsl::memory_monitor::event(key);
    auto da = sdsl::cache_file_name<std::vector<int>>(
        t_config.keys[conf::kDA].get<std::string>(), t_config);
    TSLP slp;
    construct(slp, t_config, da);
  }
};

// Bare-diff — base DifferentialSLP under kDSLPNS, at the fixed kDiffBlockSize.
template <typename TSLP, typename TRoots, typename TSpanSums, typename TSamples,
          typename TSampleRootsPos, typename TBV>
struct SlpNsTraits<DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>> {
  using TDiff = DifferentialSLP<TSLP, TRoots, TSpanSums, TSamples, TSampleRootsPos, TBV>;
  static constexpr std::string_view kKey = conf::kDSLPNS;
  static void build(Config& t_config) {
    const auto key = t_config.keys[conf::kDSLPNS].get<std::string>();
    if (sdsl::cache_file_exists<TDiff>(key, t_config)) return;
    auto event = sdsl::memory_monitor::event(key);
    TDiff dslp;
    construct(dslp, t_config, kDiffBlockSize, key);
  }
};

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TSLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>>
class DocListIdxSLP : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using StorageBase = DocListIndexExtStorage<TStorage, TAlphabet>;
  using StorageBaseImpl = IndexBaseWithExternalStorage<TStorage, TAlphabet::int_width>;
  using TPattern = typename TAlphabet::string_type;
  using TDocId = std::size_t;
  using size_type = std::size_t;

  explicit DocListIdxSLP(const TStorage& t_storage) : StorageBase(t_storage), count_idx_(t_storage) {}

  DocListIdxSLP(const TStorage& t_storage, const TCountIdx& t_count_idx)
      : StorageBase(t_storage), count_idx_(t_count_idx) {}

  DocListIdxSLP() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = count_idx_.Count(t_pattern);
    if (sp >= ep) return;

    std::vector<TDocId> docs;
    docs.reserve(ep - sp);
    auto add = [&docs](auto v) { docs.emplace_back(static_cast<TDocId>(v)); };

    // Plain SLP → span-cover expansion; differential SLP → range decompression.
    SlpNsExpandRange(*slp_, sp, ep, add);

    std::sort(docs.begin(), docs.end());
    docs.erase(std::unique(docs.begin(), docs.end()), docs.end());

    for (auto d : docs) t_report(d);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return count_idx_.serialize(out, child, "count_idx") + slp_->serialize(out, child, "slp");
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    if (slp_) append(r, "slp", sdsl::size_in_bytes(*slp_));
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    return r;
  }

  const TCountIdx& count_idx() const { return count_idx_; }

 protected:
  void loadInner(typename StorageBaseImpl::TSource& t_source, const JSON& t_keys) override {
    using namespace conf;

    std::visit([this](auto&& tt_source) { count_idx_.load(tt_source.get()); }, t_source);

    slp_ = this->template loadItemPtr<TSLP>(
        t_keys[SlpNsTraits<TSLP>::kKey].template get<std::string>(), t_source, true);
  }

  TCountIdx count_idx_;
  const TSLP* slp_ = nullptr;
};

//~~~~~~~


// Build a bare grammar SLP from a DA file via RePair. Mirrors the RePair
// step in `construct(grammar::CombinedSLP&, ...)` but without the sampled-
// tree wrapper. Stores under `t_keys[kSLPNS]` (type-hashed, so distinct
// from any GCDA-side SLP cache file).
template <typename TVarsContainer, typename TLengthsContainer>
void construct(grammar::SLP<TVarsContainer, TLengthsContainer>& t_slp,
               Config& t_config,
               const std::string& t_datafile);


// Top-level construct for the non-sampled SLP index. Same boilerplate as
// `construct(DocListIdxGCDA&, Config&)` for kText/kSA/kDocEnds/kDA, but
// builds a bare SLP under kSLPNS instead of the GCDA-sampled SLP under
// kGCDA::kSLP.
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TSLP>
void construct(DocListIdxSLP<TStorage, TAlphabet, TCountIdx, TSLP>& t_index, Config& t_config) {
  using namespace conf;

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

  // Match GCDA's cache-existence check: ConstructDocArray writes DA with a
  // type hash, so a no-hash check would never see it.
  if (const auto key = t_config.keys[kDA].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray(t_config);
  }

  // Build the SLP cache — plain bare grammar::SLP under kSLPNS, or the base
  // differential SLP (bare-diff) under kDSLPNS — dispatched by TSLP.
  SlpNsTraits<TSLP>::build(t_config);

  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  t_index.load(t_config);
}

//~~~~~~~


template <typename TVarsContainer, typename TLengthsContainer>
void construct(grammar::SLP<TVarsContainer, TLengthsContainer>& t_slp,
               Config& t_config,
               const std::string& t_datafile) {
  using namespace conf;

  // Run RePair on the DA file (only if the .R/.C outputs aren't already
  // present from an earlier sibling build — the file pair is shared
  // across all SLP-using variants and is independent of TSLP container types).
  if (!std::filesystem::exists(t_datafile + ".R") && repair::kAvailable) {
    const auto filename = std::filesystem::path(t_datafile).filename().string();
    auto event = sdsl::memory_monitor::event("RePair-" + filename);
    RunRePair(t_datafile, t_config.repair);
    t_config.file_map[filename + ".R"] = t_datafile + ".R";
    t_config.file_map[filename + ".C"] = t_datafile + ".C";
  }
  CheckRePairGrammar(t_datafile);

  auto key_slp = t_config.keys[kSLPNS].get<std::string>();
  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::SLP<TVarsContainer, TLengthsContainer>>(key_slp, t_config));

  // Build into a default-typed SLP first (push_back-friendly std::vector<uint32_t>
  // containers), then convert via the SLP template copy constructor — sdsl::dac_vector
  // and sdsl::vlc_vector have no push_back, so RePairReader cannot fill them
  // directly. The conversion uses grammar::Construct overloads (grammar/utility.h)
  // which dispatch via SFINAE: dac_vector / vlc_vector use their templated Container
  // constructor (which bit-packs internally); sdsl::int_vector falls through to the
  // resize-and-std::copy branch (which leaves the default 64-bit width). Pass a
  // bit-compress action so int_vector targets get tight per-element width; the
  // action no-ops on dac_vector / vlc_vector via `if constexpr` detection.
  grammar::SLP<> tmp_slp;
  {
    grammar::RePairReader<true> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(tmp_slp);
    re_pair_reader.Read(t_datafile, slp_wrapper);
  }
  auto bc = []<typename TVec>(TVec& v) {
    // Probe for bit_resize/width members directly — `requires { bit_compress(vv); }`
    // would unhelpfully pass for any T (since bit_compress is an unconstrained
    // template) and only fail at body instantiation. int_vector has both members;
    // std::vector and dac_vector / vlc_vector have neither.
    if constexpr (requires(TVec& vv) {
                    vv.bit_resize(std::size_t{});
                    vv.width(uint8_t{});
                  }) {
      sdsl::util::bit_compress(v);
    }
  };
  t_slp = grammar::SLP<TVarsContainer, TLengthsContainer>(tmp_slp, bc, bc);

  sdsl::store_to_cache(t_slp, key_slp, t_config, true);
}

}  // namespace dret
