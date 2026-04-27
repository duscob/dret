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
#include <sdsl/memory_management.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "sr-index/sr_idx_generic.h"
#include "sr-index/sr_index.h"

#include "construct_base.h"
#include "doc_list_index.h"
#include "index_base.h"

namespace dret {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::SrIdxGeneric<sri::SrIndexValidArea<TStorage, TAlphabet>, 16>,
          typename TSLP = grammar::SLP<>>
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

    std::vector<typename TSLP::VariableType> cover;
    grammar::ComputeSpanCover(*slp_, sp, ep, std::back_inserter(cover));

    std::vector<TDocId> docs;
    docs.reserve(ep - sp);
    auto add = [&docs](auto v) { docs.emplace_back(static_cast<TDocId>(v)); };

    for (auto var : cover) {
      auto length = slp_->SpanLength(var);
      grammar::ExpandSLPForward(slp_->GetRules(), slp_->Sigma(), var, length, add);
    }

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

    slp_ = this->template loadItemPtr<TSLP>(t_keys[kSLPNS].get<std::string>(), t_source, true);
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

  if (const auto key = t_config.keys[kSLPNS].get<std::string>(); !sdsl::cache_file_exists<TSLP>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    auto filepath_da = sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config);
    TSLP slp;
    construct(slp, t_config, filepath_da);
  }

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
  // present from an earlier GCDA-side build — the file pair is shared
  // across all SLP-using variants).
  if (!std::filesystem::exists(t_datafile + ".R") && REPAIR_EXE) {
    const auto filename = std::filesystem::path(t_datafile).filename().string();
    auto event = sdsl::memory_monitor::event("RePair-" + filename);
    std::string cmd = REPAIR_EXE + (" " + t_datafile);
    std::system(cmd.c_str());
    t_config.file_map[filename + ".R"] = t_datafile + ".R";
    t_config.file_map[filename + ".C"] = t_datafile + ".C";
  }

  auto key_slp = t_config.keys[kSLPNS].get<std::string>();
  auto event = sdsl::memory_monitor::event(
      sdsl::cache_file_name<grammar::SLP<TVarsContainer, TLengthsContainer>>(key_slp, t_config));

  {
    grammar::RePairReader<true> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(t_slp);
    re_pair_reader.Read(t_datafile, slp_wrapper);
  }

  sdsl::store_to_cache(t_slp, key_slp, t_config, true);
}

}  // namespace dret
