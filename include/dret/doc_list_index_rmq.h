//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// RMQ-based document listing index (SADA / ILCP / CILCP variants), mirroring the
// frequency-index structure in doc_freq_index_rmq.h but trimmed for pure listing:
// only the LEFTMOST pass is needed to visit each distinct doc in [sp, ep] once.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/sd_vector.hpp>

#include "sr-index/sr_idx_generic.h"
#include "sr-index/sr_index.h"

#include "config.h"
#include "construct_base.h"
#include "doc_index_rmq_common.h"
#include "doc_list_index.h"
#include "index_base.h"
#include "size_report.h"

namespace dret {
namespace rmq {

//~~~~~~~
// Internals


namespace internal {

// Number of documents = number of '\1' delimiters in the cached text (which is the
// post-`ConstructText` form: documents separated by `\1`, followed by a final `\0`).
inline std::size_t ReadNDoc(const Config& t_config) {
  sdsl::int_vector_buffer<8> text_buf(
      sdsl::cache_file_name(t_config.keys[conf::kText].get<std::string>(), t_config));
  std::size_t n = 0;
  for (std::size_t i = 0; i < text_buf.size(); ++i) {
    if (text_buf[i] == 1) ++n;
  }
  return n;
}

// Build the backward-ILCP array: for each SA position i with document d and local
// doc-suffix index j, ILCP[i] = LCP between that suffix and the previous suffix of
// the same doc (0 when j == 0). Positions where DA[i] is the sentinel terminator
// doc (DA[i] >= n_doc) are assigned 0; they never appear in a query's [sp, ep].
template <uint8_t t_width>
sdsl::int_vector<> ComputeIlcpBackward(Config& t_config, std::size_t t_n_doc) {
  using namespace dret::conf;

  // Per-doc LCP arrays.
  std::vector<sdsl::int_vector<>> doc_lcps;
  doc_lcps.reserve(t_n_doc);

  constexpr uint8_t kInternalDocDelim = 1;

  {
    sdsl::int_vector_buffer<t_width> text_buf(
        sdsl::cache_file_name(t_config.keys[kText].get<std::string>(), t_config));

    sdsl::cache_config doc_config(false, t_config.dir, t_config.id + "-doc");

    std::size_t pos = 0;
    while (pos < text_buf.size() && doc_lcps.size() < t_n_doc) {
      sdsl::int_vector<t_width> text_doc;
      {
        std::vector<typename sdsl::int_vector_buffer<t_width>::value_type> buf;
        for (; pos < text_buf.size(); ++pos) {
          auto ch = text_buf[pos];
          if (ch == kInternalDocDelim) break;
          buf.push_back(ch);
        }
        ++pos;  // skip delimiter
        text_doc.resize(buf.size() + 1);
        std::copy(buf.begin(), buf.end(), text_doc.begin());
        text_doc[text_doc.size() - 1] = 0;
      }

      sdsl::store_to_cache(text_doc, sdsl::conf::KEY_TEXT, doc_config);
      sdsl::construct_sa<t_width>(doc_config);
      sdsl::construct_isa(doc_config);
      sdsl::construct_lcp_kasai<t_width>(doc_config);

      sdsl::int_vector<> lcp;
      sdsl::load_from_cache(lcp, sdsl::conf::KEY_LCP, doc_config);
      doc_lcps.emplace_back(std::move(lcp));

      for (const auto& item : doc_config.file_map) {
        std::remove(item.second.c_str());
      }
      doc_config.file_map.clear();
    }
  }

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  // Pre-allocate at a safe width (64 bits/element) so writes never truncate; we
  // bit-compress at the end. Using the default width of 1 here would truncate values
  // to 0/1 and (because int_vector<>::width() retroactively reinterprets the same
  // bits) leave size() inconsistent with the element count.
  sdsl::int_vector<> ilcp(da.size(), 0, 64);
  std::vector<std::size_t> doc_suffix_indices(doc_lcps.size(), 0);

  for (std::size_t i = 0; i < da.size(); ++i) {
    std::size_t doc = da[i];
    if (doc >= doc_lcps.size()) {
      ilcp[i] = 0;
      continue;
    }
    auto& idx = doc_suffix_indices[doc];
    const auto& doc_lcp = doc_lcps[doc];
    ilcp[i] = doc_lcp[idx];
    ++idx;
  }

  sdsl::util::bit_compress(ilcp);
  return ilcp;
}

}  // namespace internal

//~~~~~~~
// SadaCore


template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>>
class SadaCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit SadaCore(const TStorage& t_storage) : Base(t_storage) {}
  SadaCore() = default;

  template <typename TReport>
  void findDocs(std::size_t t_sp, std::size_t t_ep, const TReport& t_report) const {
    if (t_sp > t_ep || !rmq_ || !da_) return;

    sdsl::bit_vector marked(n_doc_, 0);

    auto set_first = [](std::size_t sp, std::size_t ep, auto& stack, OccurrenceSide) {
      stack.emplace(sp, ep);
    };
    auto get_doc = [this](std::size_t i, OccurrenceSide, std::size_t, std::size_t) {
      return static_cast<std::size_t>((*da_)[i]);
    };
    auto is_reported = [&marked](std::size_t, std::size_t doc, OccurrenceSide) {
      return doc >= marked.size() || marked[doc];
    };
    auto report = [&marked, &t_report](
                      std::size_t, std::size_t doc, OccurrenceSide, std::size_t, std::size_t) {
      marked[doc] = 1;
      t_report(doc);
    };

    GetExtremeOccurrencesRMQ<OccurrenceSide::LEFTMOST>(
        t_sp, t_ep, set_first, *rmq_, get_doc, is_reported, report);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += rmq_ ? sdsl::serialize(*rmq_, out, child, "rmq") : sdsl::serialize_empty_object<TRMQ>(out, child, "rmq");
    written += da_ ? sdsl::serialize(*da_, out, child, "da")
                   : sdsl::serialize_empty_object<sdsl::int_vector<>>(out, child, "da");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rmq_) append(r, "rmq", sdsl::size_in_bytes(*rmq_));
    if (da_) append(r, "da", sdsl::size_in_bytes(*da_));
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    key_rmq_ = t_keys[kSADA][kRmq].get<std::string>();
    key_da_ = t_keys[kDA].get<std::string>();

    rmq_ = this->template loadItemPtr<TRMQ>(key_rmq_, t_source, true);
    da_ = this->template loadItemPtr<sdsl::int_vector<>>(key_da_, t_source, true);

    // n_doc bound: max(DA) + 1 — over-counts by one when the terminator suffix is
    // present, but `is_reported` filters out-of-range docs anyway.
    n_doc_ = 0;
    for (std::size_t i = 0; i < da_->size(); ++i) {
      auto v = static_cast<std::size_t>((*da_)[i]);
      if (v + 1 > n_doc_) n_doc_ = v + 1;
    }
  }

  std::string key_rmq_;
  std::string key_da_;

  const TRMQ* rmq_ = nullptr;
  const sdsl::int_vector<>* da_ = nullptr;
  std::size_t n_doc_ = 0;
};

//~~~~~~~
// IlcpCore / CilcpCore (share the same loaded representation; differ in the
// run-compression rule used at construction and in the fan-out behavior at query).


enum class IlcpVariant { ILCP, CILCP };

template <IlcpVariant kVariant,
          typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds>
class IlcpLikeCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit IlcpLikeCore(const TStorage& t_storage) : Base(t_storage) {}
  IlcpLikeCore() = default;

  template <typename TReport>
  void findDocs(std::size_t t_sp, std::size_t t_ep, const TReport& t_report) const {
    if (t_sp > t_ep || !rmq_ || !da_ || !run_heads_ || !rank_ || !select_) return;

    sdsl::bit_vector marked(n_doc_, 0);
    const auto& rank = *rank_;
    const auto& select = *select_;

    std::size_t run_sp = rank(t_sp + 1) - 1;
    std::size_t run_ep = rank(t_ep + 1) - 1;

    auto set_first = [run_sp, run_ep](std::size_t, std::size_t, auto& stack, OccurrenceSide) {
      stack.emplace(run_sp, run_ep);
    };
    auto get_doc =
        [this, &select](std::size_t i, OccurrenceSide, std::size_t sp_orig, std::size_t) {
          std::size_t head = static_cast<std::size_t>(select(i + 1));
          std::size_t pos = std::max(sp_orig, head);
          return static_cast<std::size_t>((*da_)[pos]);
        };
    auto is_reported = [&marked](std::size_t, std::size_t doc, OccurrenceSide) {
      return doc >= marked.size() || marked[doc];
    };
    auto report = [this, &marked, &t_report, &select](std::size_t i,
                                                       std::size_t doc,
                                                       OccurrenceSide,
                                                       std::size_t sp_orig,
                                                       std::size_t ep_orig) {
      marked[doc] = 1;
      t_report(doc);

      // Fan out the rest of the run in [head+1, tail] ∩ [sp, ep]. For the last
      // run, select(i+2) is undefined; fall back to da.size() as the next-head.
      std::size_t head = static_cast<std::size_t>(select(i + 1));
      std::size_t run_start = std::max(sp_orig, head);
      std::size_t next_head = (i + 1 < n_runs_) ? static_cast<std::size_t>(select(i + 2)) : da_->size();
      std::size_t run_end = std::min(ep_orig, next_head - 1);
      for (std::size_t p = run_start + 1; p <= run_end; ++p) {
        auto d = static_cast<std::size_t>((*da_)[p]);
        if constexpr (kVariant == IlcpVariant::CILCP) {
          if (d == doc) break;  // CILCP: stop when we loop back to the run's anchor doc
        }
        if (d < marked.size() && !marked[d]) {
          marked[d] = 1;
          t_report(d);
        }
      }
    };

    GetExtremeOccurrencesRMQ<OccurrenceSide::LEFTMOST>(
        t_sp, t_ep, set_first, *rmq_, get_doc, is_reported, report);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += rmq_ ? sdsl::serialize(*rmq_, out, child, "rmq")
                    : sdsl::serialize_empty_object<TRMQ>(out, child, "rmq");
    written += run_heads_ ? sdsl::serialize(*run_heads_, out, child, "run_heads")
                          : sdsl::serialize_empty_object<TBvRunHeads>(out, child, "run_heads");
    written += rank_ ? sdsl::serialize(*rank_, out, child, "run_heads_rank")
                     : sdsl::serialize_empty_object<typename TBvRunHeads::rank_1_type>(out, child, "run_heads_rank");
    written += select_
                   ? sdsl::serialize(*select_, out, child, "run_heads_select")
                   : sdsl::serialize_empty_object<typename TBvRunHeads::select_1_type>(out, child, "run_heads_select");
    written += da_ ? sdsl::serialize(*da_, out, child, "da")
                   : sdsl::serialize_empty_object<sdsl::int_vector<>>(out, child, "da");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rmq_) append(r, "rmq", sdsl::size_in_bytes(*rmq_));
    if (run_heads_) append(r, "run_heads", sdsl::size_in_bytes(*run_heads_));
    if (da_) append(r, "da", sdsl::size_in_bytes(*da_));
    return r;
  }

 protected:
  using typename Base::TSource;

  static std::string_view TopKey() {
    return kVariant == IlcpVariant::ILCP ? conf::kILCP : conf::kCILCP;
  }

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    std::string top(TopKey());
    key_rmq_ = t_keys[top][kRmq].get<std::string>();
    key_run_heads_ = t_keys[top][kRunHeads].get<std::string>();
    key_da_ = t_keys[kDA].get<std::string>();
    auto key_doc_end = t_keys[kDocEnds].get<std::string>();

    rmq_ = this->template loadItemPtr<TRMQ>(key_rmq_, t_source, true);
    run_heads_ = this->template loadItemPtr<TBvRunHeads>(key_run_heads_, t_source, true);
    rank_ = &this->template loadBVRank<TBvRunHeads>(key_run_heads_, t_source, true).get();
    select_ = &this->template loadBVSelect<TBvRunHeads>(key_run_heads_, t_source, true).get();
    da_ = this->template loadItemPtr<sdsl::int_vector<>>(key_da_, t_source, true);

    n_runs_ = (*rank_)(run_heads_->size());

    // Derive n_doc from the DA max (avoids rank_1(size) edge-case on sd_vector).
    n_doc_ = 0;
    for (std::size_t i = 0; i < da_->size(); ++i) {
      auto v = static_cast<std::size_t>((*da_)[i]);
      if (v + 1 > n_doc_) n_doc_ = v + 1;
    }
  }

  std::string key_rmq_;
  std::string key_run_heads_;
  std::string key_da_;

  const TRMQ* rmq_ = nullptr;
  const TBvRunHeads* run_heads_ = nullptr;
  const typename TBvRunHeads::rank_1_type* rank_ = nullptr;
  const typename TBvRunHeads::select_1_type* select_ = nullptr;
  const sdsl::int_vector<>* da_ = nullptr;
  std::size_t n_doc_ = 0;
  std::size_t n_runs_ = 0;
};

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TBvRunHeads = sdsl::bit_vector,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>>
using IlcpCore = IlcpLikeCore<IlcpVariant::ILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds>;

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TBvRunHeads = sdsl::bit_vector,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>>
using CilcpCore = IlcpLikeCore<IlcpVariant::CILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds>;

//~~~~~~~
// DocListIdxRMQ


template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::SrIdxGeneric<sri::SrIndexValidArea<TStorage, TAlphabet>, 16>,
          typename TCore = SadaCore<TStorage, TAlphabet::int_width>>
class DocListIdxRMQ : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using Base = DocListIndexExtStorage<TStorage, TAlphabet>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  explicit DocListIdxRMQ(const TStorage& t_storage) : Base(t_storage), count_idx_(t_storage), core_(t_storage) {}

  DocListIdxRMQ(const TStorage& t_storage, const TCountIdx& t_count_idx)
      : Base(t_storage), count_idx_(t_count_idx), core_(t_storage) {}

  DocListIdxRMQ() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = count_idx_.Count(t_pattern);
    if (sp >= ep) return;          // sri::Count uses half-open [sp, ep)
    core_.findDocs(sp, ep - 1, t_report);  // convert to closed [sp, ep] for RMQ
  }

  void load(Config t_config) override {
    count_idx_.load(t_config);
    core_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    count_idx_.load(in);
    core_.load(in, t_keys);
  }

  using Base::load;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    return count_idx_.serialize(out, child, "count_idx") + core_.serialize(out, child, "core");
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    auto core_report = core_.GetSizeReport();
    for (const auto& [k, v] : core_report) append(r, k, v);
    return r;
  }

  const TCountIdx& count_idx() const {
    return count_idx_;
  }

 protected:
  TCountIdx count_idx_;
  TCore core_;
};

//~~~~~~~
// Construction


namespace internal {

// Ensure Text/SA/DocEnds/DA are present in the cache. Each call is a no-op when
// the respective cache entry already exists.
template <uint8_t t_width, typename TBvDocEnds>
void EnsureBasicStructures(Config& t_config) {
  using namespace dret::conf;

  if (!cache_file_exists(t_config.keys[kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<t_width>(t_config);
  }

  if (const auto key = t_config.keys[kSA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    sdsl::construct_sa<t_width>(t_config);
  }

  if (const auto key = t_config.keys[kDocEnds].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocEnd<t_width, TBvDocEnds>(t_config);
  }

  // ConstructDocArray writes DA only with type-hash; the no-hash check would never
  // see it, causing DA to be rebuilt on every construct() call. Match the storage.
  if (const auto key = t_config.keys[kDA].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray<TBvDocEnds>(t_config);
  }

  if (const auto key = t_config.keys[kRmqNDoc].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<64>>(key, t_config)) {
    sdsl::int_vector<64> n_doc_vec(1, ReadNDoc(t_config));
    sdsl::store_to_cache(n_doc_vec, key, t_config, true);
  }
}

// Store run_heads bitvector (with rank/select) and the RMQ built from run values.
template <typename TBvRunHeads, typename TRMQ>
void StoreRunHeadsAndRMQ(Config& t_config,
                         const std::string& t_key_run_heads,
                         const std::string& t_key_rmq,
                         sdsl::bit_vector&& t_run_heads_bv,
                         std::vector<std::size_t>& t_run_values) {
  TBvRunHeads run_heads(t_run_heads_bv);
  sdsl::store_to_cache(run_heads, t_key_run_heads, t_config, true);

  typename TBvRunHeads::rank_1_type rank(&run_heads);
  sdsl::store_to_cache(rank, t_key_run_heads, t_config, true);

  typename TBvRunHeads::select_1_type select(&run_heads);
  sdsl::store_to_cache(select, t_key_run_heads, t_config, true);

  TRMQ rmq(&t_run_values);
  sdsl::store_to_cache(rmq, t_key_rmq, t_config, true);
}

}  // namespace internal

// SADA construction: prev_doc + RMinQ.
template <typename TStorage, uint8_t t_width, typename TRMQ, typename TBvDocEnds>
void construct(SadaCore<TStorage, t_width, TRMQ, TBvDocEnds>& t_core, Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kSADA][kRmq].get<std::string>();
  if (sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)) return;

  auto event = sdsl::memory_monitor::event(key_rmq);

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  std::size_t n_doc = internal::ReadNDoc(t_config);

  sdsl::int_vector<> prev_doc(da.size(), 0, sdsl::bits::hi(da.size()) + 1);
  std::vector<std::size_t> last_occ(n_doc + 2, 0);
  for (std::size_t i = 0; i < da.size(); ++i) {
    std::size_t doc = da[i];
    if (doc >= last_occ.size()) last_occ.resize(doc + 1, 0);
    prev_doc[i] = last_occ[doc];
    last_occ[doc] = i;
  }

  TRMQ rmq(&prev_doc);
  sdsl::store_to_cache(rmq, key_rmq, t_config, true);
}

// ILCP construction: plain RLE on backward-ILCP.
template <typename TStorage, uint8_t t_width, typename TBvRunHeads, typename TRMQ, typename TBvDocEnds>
void construct(
    IlcpLikeCore<IlcpVariant::ILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds>& t_core,
    Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kILCP][kRmq].get<std::string>();
  auto key_run_heads = t_config.keys[kILCP][kRunHeads].get<std::string>();
  if (sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      && sdsl::cache_file_exists<TBvRunHeads>(key_run_heads, t_config)) {
    return;
  }

  auto event = sdsl::memory_monitor::event(key_rmq);

  std::size_t n_doc = internal::ReadNDoc(t_config);
  auto ilcp = internal::ComputeIlcpBackward<t_width>(t_config, n_doc);

  sdsl::bit_vector run_heads(ilcp.size(), 0);
  std::vector<std::size_t> run_values;
  run_values.emplace_back(ilcp[0]);
  run_heads[0] = 1;
  for (std::size_t i = 1; i < ilcp.size(); ++i) {
    if (ilcp[i - 1] != ilcp[i]) {
      run_values.emplace_back(ilcp[i]);
      run_heads[i] = 1;
    }
  }

  internal::StoreRunHeadsAndRMQ<TBvRunHeads, TRMQ>(
      t_config, key_run_heads, key_rmq, std::move(run_heads), run_values);
}

// CILCP construction: DA-aware RLE rule on backward-ILCP.
template <typename TStorage, uint8_t t_width, typename TBvRunHeads, typename TRMQ, typename TBvDocEnds>
void construct(
    IlcpLikeCore<IlcpVariant::CILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds>& t_core,
    Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kCILCP][kRmq].get<std::string>();
  auto key_run_heads = t_config.keys[kCILCP][kRunHeads].get<std::string>();
  if (sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      && sdsl::cache_file_exists<TBvRunHeads>(key_run_heads, t_config)) {
    return;
  }

  auto event = sdsl::memory_monitor::event(key_rmq);

  std::size_t n_doc = internal::ReadNDoc(t_config);
  auto ilcp = internal::ComputeIlcpBackward<t_width>(t_config, n_doc);

  sdsl::int_vector<> da;
  sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

  sdsl::bit_vector run_heads(ilcp.size(), 0);
  std::vector<std::size_t> run_values;
  run_values.emplace_back(ilcp[0]);
  run_heads[0] = 1;
  for (std::size_t i = 1; i < ilcp.size(); ++i) {
    std::size_t l = ilcp[i - 1];
    std::size_t d = da[i - 1];
    if (l == ilcp[i]) {
      while (++i < ilcp.size() && l == ilcp[i]) {
      }
    } else if (d == da[i]) {
      do {
        l = std::min<std::size_t>(l, ilcp[i]);
      } while (++i < ilcp.size() && d == da[i]);
    }

    if (i < ilcp.size()) {
      run_values[run_values.size() - 1] = l;
      run_values.emplace_back(ilcp[i]);
      run_heads[i] = 1;
    }
  }

  internal::StoreRunHeadsAndRMQ<TBvRunHeads, TRMQ>(
      t_config, key_run_heads, key_rmq, std::move(run_heads), run_values);
}

// Wiring: construct the whole DocListIdxRMQ (count_idx + core).
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TCore>
void construct(DocListIdxRMQ<TStorage, TAlphabet, TCountIdx, TCore>& t_index, Config& t_config) {
  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  TCore core(t_index.storage());
  construct(core, t_config);

  t_index.load(t_config);
}

}  // namespace rmq
}  // namespace dret
