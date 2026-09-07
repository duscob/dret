//
// Created by Dustin Cobas <dustin.cobas@gmail.com>.
//
// RMQ-based document listing index (SADA / ILCP / CILCP variants), mirroring the
// frequency-index structure in doc_freq/doc_freq_rmq.h but trimmed for pure listing:
// only the LEFTMOST pass is needed to visit each distinct doc in [sp, ep] once.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/sd_vector.hpp>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/rmq/doc_index_rmq_common.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/rmq/doc_list_rmq_scheme.h"
#include "dret/index_base.h"
#include "dret/rmq/rmq_get_doc_policies.h"
#include "dret/size_report.h"

namespace dret {
namespace rmq {

//~~~~~~~
// Internals


namespace internal {

// Number of documents = number of '\1' delimiters in the cached text (which is the
// post-`ConstructText` form: documents separated by `\1`, followed by a final `\0`).
inline std::size_t ReadNDoc(const Config& t_config) {
  sdsl::int_vector_buffer<8> text_buf(sdsl::cache_file_name(t_config.keys[conf::kText].get<std::string>(), t_config));
  std::size_t n = 0;
  for (std::size_t i = 0; i < text_buf.size(); ++i) {
    if (text_buf[i] == 1)
      ++n;
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
    sdsl::int_vector_buffer<t_width> text_buf(sdsl::cache_file_name(t_config.keys[kText].get<std::string>(), t_config));

    sdsl::cache_config doc_config(false, t_config.dir, t_config.id + "-doc");

    std::size_t pos = 0;
    while (pos < text_buf.size() && doc_lcps.size() < t_n_doc) {
      sdsl::int_vector<t_width> text_doc;
      {
        std::vector<typename sdsl::int_vector_buffer<t_width>::value_type> buf;
        for (; pos < text_buf.size(); ++pos) {
          auto ch = text_buf[pos];
          if (ch == kInternalDocDelim)
            break;
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

// Load the cached backward-ILCP array, or compute it once and cache it. The array
// is collection-level (independent of core, get-doc, run encoding), so all four
// ILCP-family cores (ILCP / CILCP / ILCP-S / CILCP-S) share one computation
// instead of recomputing the expensive per-document LCP up to 4x on a cold build.
template <uint8_t t_width>
sdsl::int_vector<> LoadOrComputeIlcp(Config& t_config, std::size_t t_n_doc) {
  const auto key = t_config.keys[conf::kIlcpArray].get<std::string>();
  sdsl::int_vector<> ilcp;
  if (sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    sdsl::load_from_cache(ilcp, key, t_config, true);
    return ilcp;
  }
  ilcp = ComputeIlcpBackward<t_width>(t_config, t_n_doc);
  sdsl::store_to_cache(ilcp, key, t_config, true);
  return ilcp;
}

}  // namespace internal

//~~~~~~~
// SadaLCore


template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>,
          typename TGetDoc = GetDocDA<TStorage, t_width>>
class SadaLCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit SadaLCore(const TStorage& t_storage) : Base(t_storage), get_doc_(t_storage) {}

  SadaLCore(const TStorage& t_storage, uint32_t t_block_size, float t_storing_factor)
      : Base(t_storage), get_doc_(t_storage, t_block_size, t_storing_factor) {}

  SadaLCore() = default;

  void load(Config t_config) override {
    Base::load(t_config);
    get_doc_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    Base::load(in, t_keys);
    get_doc_.load(in, t_keys);
  }

  using Base::load;

  template <typename TReport>
  void findDocs(std::size_t t_sp, std::size_t t_ep, std::size_t /*t_m*/, const TReport& t_report) const {
    if (t_sp >= t_ep || !rmq_)
      return;

    MarkedReported mr(n_doc_);
    auto report = [&mr, &t_report](std::size_t /*k*/, std::size_t d) {
      mr.mark(d);
      t_report(d);
    };

    ListDocsRMQScheme(t_sp, t_ep, *rmq_, get_doc_, mr, report);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += rmq_ ? sdsl::serialize(*rmq_, out, child, "rmq") : sdsl::serialize_empty_object<TRMQ>(out, child, "rmq");
    // Write n_doc so the istream deserialization path can read it in sequence.
    sdsl::int_vector<64> n_doc_vec(1, n_doc_);
    written += sdsl::serialize(n_doc_vec, out, child, "n_doc");
    written += get_doc_.serialize(out, child, "get_doc");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rmq_)
      append(r, "rmq", sdsl::size_in_bytes(*rmq_));
    auto get_doc_report = get_doc_.GetSizeReport();
    for (const auto& [k, v] : get_doc_report)
      append(r, k, v);
    return r;
  }

  TGetDoc& get_doc_policy() {
    return get_doc_;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    key_rmq_ = t_keys[kSADA][kRmq].get<std::string>();
    rmq_ = this->template loadItemPtr<TRMQ>(key_rmq_, t_source, true);

    // Config path: loads from the kRmqNDoc cache entry (key-based lookup).
    // Istream path: reads the next int_vector<64> sequentially from the stream,
    // matching the n_doc_vec written by serialize() above.
    const auto key_n_doc = t_keys[kRmqNDoc].get<std::string>();
    const auto* n_doc_vec = this->template loadItemPtr<sdsl::int_vector<64>>(key_n_doc, t_source, true);
    n_doc_ = static_cast<std::size_t>((*n_doc_vec)[0]);
  }

  std::string key_rmq_;

  const TRMQ* rmq_ = nullptr;
  std::size_t n_doc_ = 0;
  TGetDoc get_doc_;
};

//~~~~~~~
// SadaCore — Sadakane-style depth-based stop variant of SadaLCore.
//
// Same RMQ over prev_doc as the original SADA core; the difference is at
// query time only. SadaCore persists the prev_doc array (existing SadaLCore
// discards it after RMQ construction) so the recursion stop can read it as
// `prev_doc[k] >= bp_subrange` — Sadakane's canonical predicate. Reuses the
// `kSADA / kRmq` cache file; adds `kSADA / kPrevDoc` for the value array.


template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>,
          typename TGetDoc = GetDocDA<TStorage, t_width>,
          typename TPrevDoc = sdsl::int_vector<>>
class SadaCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit SadaCore(const TStorage& t_storage) : Base(t_storage), get_doc_(t_storage) {}

  SadaCore(const TStorage& t_storage, uint32_t t_block_size, float t_storing_factor)
      : Base(t_storage), get_doc_(t_storage, t_block_size, t_storing_factor) {}

  SadaCore() = default;

  void load(Config t_config) override {
    Base::load(t_config);
    get_doc_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    Base::load(in, t_keys);
    get_doc_.load(in, t_keys);
  }

  using Base::load;

  template <typename TReport>
  void findDocs(std::size_t t_sp, std::size_t t_ep, std::size_t /*t_m*/, const TReport& t_report) const {
    if (t_sp >= t_ep || !rmq_ || !prev_doc_)
      return;

    MarkedReported mr(n_doc_);
    auto report = [&mr, &t_report](std::size_t /*k*/, std::size_t d) {
      mr.mark(d);
      t_report(d);
    };
    auto stop_pred = [this](std::size_t k, std::size_t bp) {
      return static_cast<std::size_t>((*prev_doc_)[k]) >= bp;
    };

    ListDocsRMQSchemeDepth(t_sp, t_ep, *rmq_, get_doc_, stop_pred, mr, report);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += rmq_ ? sdsl::serialize(*rmq_, out, child, "rmq") : sdsl::serialize_empty_object<TRMQ>(out, child, "rmq");
    written += prev_doc_
                   ? sdsl::serialize(*prev_doc_, out, child, "prev_doc")
                   : sdsl::serialize_empty_object<TPrevDoc>(out, child, "prev_doc");
    sdsl::int_vector<64> n_doc_vec(1, n_doc_);
    written += sdsl::serialize(n_doc_vec, out, child, "n_doc");
    written += get_doc_.serialize(out, child, "get_doc");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rmq_)
      append(r, "rmq", sdsl::size_in_bytes(*rmq_));
    if (prev_doc_)
      append(r, "prev_doc", sdsl::size_in_bytes(*prev_doc_));
    auto get_doc_report = get_doc_.GetSizeReport();
    for (const auto& [k, v] : get_doc_report)
      append(r, k, v);
    return r;
  }

  TGetDoc& get_doc_policy() {
    return get_doc_;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    key_rmq_ = t_keys[kSADA][kRmq].get<std::string>();
    rmq_ = this->template loadItemPtr<TRMQ>(key_rmq_, t_source, true);

    key_prev_doc_ = t_keys[kSADA][kPrevDoc].get<std::string>();
    prev_doc_ = this->template loadItemPtr<TPrevDoc>(key_prev_doc_, t_source, true);

    const auto key_n_doc = t_keys[kRmqNDoc].get<std::string>();
    const auto* n_doc_vec = this->template loadItemPtr<sdsl::int_vector<64>>(key_n_doc, t_source, true);
    n_doc_ = static_cast<std::size_t>((*n_doc_vec)[0]);
  }

  std::string key_rmq_;
  std::string key_prev_doc_;

  const TRMQ* rmq_ = nullptr;
  const TPrevDoc* prev_doc_ = nullptr;
  std::size_t n_doc_ = 0;
  TGetDoc get_doc_;
};

//~~~~~~~
// IlcpLCore / CilcpLCore (share the same loaded representation; differ in the
// run-compression rule used at construction and in the fan-out behavior at query).


enum class IlcpLeanVariant { ILCP_L, CILCP_L };

template <IlcpLeanVariant kVariant,
          typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds,
          typename TGetDoc = GetDocDA<TStorage, t_width>>
class IlcpLikeLeanCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit IlcpLikeLeanCore(const TStorage& t_storage) : Base(t_storage), get_doc_(t_storage) {}

  IlcpLikeLeanCore(const TStorage& t_storage, uint32_t t_block_size, float t_storing_factor)
      : Base(t_storage), get_doc_(t_storage, t_block_size, t_storing_factor) {}

  IlcpLikeLeanCore() = default;

  void load(Config t_config) override {
    Base::load(t_config);
    get_doc_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    Base::load(in, t_keys);
    get_doc_.load(in, t_keys);
  }

  using Base::load;

  template <typename TReport>
  void findDocs(std::size_t t_sp, std::size_t t_ep, std::size_t /*t_m*/, const TReport& t_report) const {
    if (t_sp >= t_ep || !rmq_ || !run_heads_ || !rank_ || !select_)
      return;

    // IlcpState holds the original SA-space [sp, ep) so that get_doc and the
    // fan-out report can recover per-run SA boundaries after PreprocessILCP
    // converts the working range to run-space.
    struct IlcpState {
      const typename TBvRunHeads::rank_1_type& rank_ref;
      std::size_t sp_orig = 0;
      std::size_t ep_orig = 0;

      const auto& rank() const {
        return rank_ref;
      }

      void setInitialRange(std::size_t sp, std::size_t ep) {
        sp_orig = sp;
        ep_orig = ep;
      }
    } state{*rank_};

    std::size_t run_sp = t_sp, run_ep = t_ep;
    PreprocessILCP<IlcpState>{state}(run_sp, run_ep);

    MarkedReported mr(n_doc_);
    const auto& select = *select_;

    auto get_doc = [this, &select, &state](std::size_t i) {
      const std::size_t head = static_cast<std::size_t>(select(i + 1));
      return get_doc_(std::max(state.sp_orig, head));
    };

    auto report = [this, &mr, &t_report, &select, &state](std::size_t i, std::size_t doc) {
      auto report_dedup = [&mr, &t_report](auto d) {
        const auto dd = static_cast<std::size_t>(d);
        if (!mr(0, dd)) {
          mr.mark(dd);
          t_report(dd);
        }
      };
      // The head doc is de-duplicated like any other: CILCP reaches report()
      // even for an already-reported head, because the fan-out below may still
      // hold the first occurrence of some other document.
      report_dedup(doc);

      // Fan out remaining positions in this run within [sp_orig, ep_orig). For
      // the last run select(i+2) is undefined; run_heads_->size() is the safe sentinel.
      const std::size_t head = static_cast<std::size_t>(select(i + 1));
      const std::size_t run_start = std::max(state.sp_orig, head);
      const std::size_t next_head = (i + 1 < n_runs_) ? static_cast<std::size_t>(select(i + 2)) : run_heads_->size();
      const std::size_t run_end = std::min(state.ep_orig - 1, next_head - 1);
      // Both ILCP and CILCP issue a single range-based SLP/DSLP descent so
      // grammar-compressed GetDoc backends pay O(|cover|*height + L) per run
      // instead of L * O(height) for L per-position descents.
      //
      // CILCP adds a peek on the first leftover position. It is sound because a
      // CILCP run is either (a) a same-doc run, where nothing past the head is
      // new, or (b) an ilcp-constant run, where the peek fires only if the
      // run's first two in-range positions share a document -- so the second is
      // not that document's first occurrence, its ILCP value is >= m, and since
      // the run is ilcp-constant EVERY position in it is a repeat occurrence
      // that some other run reports. (The converse case, first two positions in
      // different documents, does not fire the peek and fans out normally.)
      const std::size_t b = run_start + 1;
      const std::size_t e = run_end + 1;
      if constexpr (kVariant == IlcpLeanVariant::CILCP_L) {
        if (b < e && get_doc_(b) != doc)
          get_doc_(b, e, report_dedup);
      } else {
        if (b < e)
          get_doc_(b, e, report_dedup);
      }
    };

    // The marker-based stop is sound for ILCP-L everywhere: its runs are
    // value-uniform, so a run's stored value equals the value at each of its
    // positions and cannot be "borrowed" from outside a query range.
    //
    // CILCP-L merges by document, so a run's stored value is a MINIMUM. For a
    // run lying entirely inside [sp, ep) that minimum is still the minimum over
    // the run's in-range positions, and value < m still implies the run holds a
    // first occurrence -- which is what makes the marker test sound there. For
    // the at most two runs that straddle a range boundary it does not: their
    // minimum may be attained outside [sp, ep), so such a run can report a
    // repeat occurrence early and suppress a later run holding a new document.
    //
    // So: recurse with the marker stop over the INTERIOR runs, and fan the
    // boundary runs out separately.
    //
    // Once the recursion is restricted that way, kAlwaysReport is no longer
    // needed here and both variants leave it at its default. On an interior run
    // an already-reported head document forces the run's stored value to be at
    // least m -- for a merged same-document run because every position after
    // the head is an adjacent occurrence of that document and so has ILCP >= m,
    // and for an unmerged run by value-uniformity. Either way the run holds no
    // first occurrence, so its fan-out cannot yield a new document and skipping
    // it is free. The boundary runs never reach this code path; they are fanned
    // out unconditionally below.
    if constexpr (kVariant != IlcpLeanVariant::CILCP_L) {
      ListDocsRMQScheme(run_sp, run_ep, *rmq_, get_doc, mr, report);
    } else {
      const std::size_t first_head = static_cast<std::size_t>(select(run_sp + 1));
      const std::size_t last = run_ep - 1;
      const std::size_t after_last =
          (last + 1 < n_runs_) ? static_cast<std::size_t>(select(last + 2)) : run_heads_->size();
      const bool clip_lo = first_head < state.sp_orig;        // run starts before sp
      const bool clip_hi = (after_last - 1) > (state.ep_orig - 1);  // run ends after ep-1

      const std::size_t lo = run_sp + (clip_lo ? 1 : 0);
      const std::size_t hi = run_ep - (clip_hi ? 1 : 0);
      if (lo < hi)
        ListDocsRMQScheme(lo, hi, *rmq_, get_doc, mr, report);

      // AFTER the interior recursion, never before: a boundary run fanned out
      // first can mark a document whose first occurrence lies in an interior
      // run, which would then be stopped on spuriously and its other documents
      // lost. (Measured: doing this first gives wrong answers.)
      if (clip_lo)
        report(run_sp, get_doc(run_sp));
      if (clip_hi && !(clip_lo && run_sp == last))
        report(last, get_doc(last));
    }
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += rmq_ ? sdsl::serialize(*rmq_, out, child, "rmq") : sdsl::serialize_empty_object<TRMQ>(out, child, "rmq");
    written += run_heads_ ? sdsl::serialize(*run_heads_, out, child, "run_heads")
                          : sdsl::serialize_empty_object<TBvRunHeads>(out, child, "run_heads");
    written += rank_ ? sdsl::serialize(*rank_, out, child, "run_heads_rank")
                     : sdsl::serialize_empty_object<typename TBvRunHeads::rank_1_type>(out, child, "run_heads_rank");
    written += select_
                   ? sdsl::serialize(*select_, out, child, "run_heads_select")
                   : sdsl::serialize_empty_object<typename TBvRunHeads::select_1_type>(out, child, "run_heads_select");
    // Write n_doc so the istream deserialization path can read it in sequence
    // (n_runs is always re-derived from rank after loading; no need to persist it).
    sdsl::int_vector<64> n_doc_vec(1, n_doc_);
    written += sdsl::serialize(n_doc_vec, out, child, "n_doc");
    written += get_doc_.serialize(out, child, "get_doc");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rmq_)
      append(r, "rmq", sdsl::size_in_bytes(*rmq_));
    if (run_heads_)
      append(r, "run_heads", sdsl::size_in_bytes(*run_heads_));
    auto get_doc_report = get_doc_.GetSizeReport();
    for (const auto& [k, v] : get_doc_report)
      append(r, k, v);
    return r;
  }

  TGetDoc& get_doc_policy() {
    return get_doc_;
  }

 protected:
  using typename Base::TSource;

  static std::string_view TopKey() {
    return kVariant == IlcpLeanVariant::ILCP_L ? conf::kILCP : conf::kCILCP;
  }

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    std::string top(TopKey());
    key_rmq_ = t_keys[top][kRmq].get<std::string>();
    key_run_heads_ = t_keys[top][kRunHeads].get<std::string>();

    rmq_ = this->template loadItemPtr<TRMQ>(key_rmq_, t_source, true);
    run_heads_ = this->template loadItemPtr<TBvRunHeads>(key_run_heads_, t_source, true);
    rank_ = &this->template loadBVRank<TBvRunHeads>(key_run_heads_, t_source, true).get();
    select_ = &this->template loadBVSelect<TBvRunHeads>(key_run_heads_, t_source, true).get();

    // n_runs is always derivable from rank after rank is loaded.
    n_runs_ = (*rank_)(run_heads_->size());

    // Config path: loads from kRmqNDoc cache; istream path: reads from stream
    // (matching the n_doc_vec written by serialize() between select and get_doc).
    const auto key_n_doc = t_keys[kRmqNDoc].get<std::string>();
    const auto* n_doc_vec = this->template loadItemPtr<sdsl::int_vector<64>>(key_n_doc, t_source, true);
    n_doc_ = static_cast<std::size_t>((*n_doc_vec)[0]);
  }

  std::string key_rmq_;
  std::string key_run_heads_;

  const TRMQ* rmq_ = nullptr;
  const TBvRunHeads* run_heads_ = nullptr;
  const typename TBvRunHeads::rank_1_type* rank_ = nullptr;
  const typename TBvRunHeads::select_1_type* select_ = nullptr;
  std::size_t n_doc_ = 0;
  std::size_t n_runs_ = 0;
  TGetDoc get_doc_;
};

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TBvRunHeads = sdsl::sd_vector<>,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>,
          typename TGetDoc = GetDocDA<TStorage, t_width>>
using IlcpLCore = IlcpLikeLeanCore<IlcpLeanVariant::ILCP_L, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc>;

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TBvRunHeads = sdsl::sd_vector<>,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>,
          typename TGetDoc = GetDocDA<TStorage, t_width>>
using CilcpLCore = IlcpLikeLeanCore<IlcpLeanVariant::CILCP_L, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc>;

//~~~~~~~
// IlcpCore / CilcpCore — Sadakane-style depth-based recursion-stop variants.
//
// Same query shape as IlcpLikeLeanCore (run-space RMQ + Option B fan-out with the
// CILCP Rule-2 peek), but the recursion stops on the canonical
// `run_values[k] >= m` predicate from Cobas, Mäkinen, Rossi SPIRE 2020
// Lemma 2 (and Sadakane's original ILCP work, Gagie-Navarro-Puglisi 2014).
//
// ILCP and ILCP-L share the ilcp_run_heads + ilcp_rmq cache files; this core
// additionally persists the run_values array (kILCP / kRunValues) that its stop
// consults. CILCP and CILCP-L stand in the same relation: both build the
// CILCP* Definition 1 partition of Cobas, Makinen and Rossi (SPIRE 2020) via
// internal::BuildCilcpRuns and share its cache, and only this core stores the
// values (kCILCP / kRunValues).


enum class IlcpFullVariant { ILCP, CILCP };

template <IlcpFullVariant kVariantS,
          typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds,
          typename TGetDoc = GetDocDA<TStorage, t_width>,
          typename TRunValues = sdsl::dac_vector<>>
class IlcpLikeFullCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit IlcpLikeFullCore(const TStorage& t_storage) : Base(t_storage), get_doc_(t_storage) {}

  IlcpLikeFullCore(const TStorage& t_storage, uint32_t t_block_size, float t_storing_factor)
      : Base(t_storage), get_doc_(t_storage, t_block_size, t_storing_factor) {}

  IlcpLikeFullCore() = default;

  void load(Config t_config) override {
    Base::load(t_config);
    get_doc_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    Base::load(in, t_keys);
    get_doc_.load(in, t_keys);
  }

  using Base::load;

  template <typename TReport>
  void findDocs(std::size_t t_sp, std::size_t t_ep, std::size_t t_m, const TReport& t_report) const {
    if (t_sp >= t_ep || !rmq_ || !run_heads_ || !rank_ || !select_ || !run_values_)
      return;

    struct IlcpState {
      const typename TBvRunHeads::rank_1_type& rank_ref;
      std::size_t sp_orig = 0;
      std::size_t ep_orig = 0;

      const auto& rank() const {
        return rank_ref;
      }

      void setInitialRange(std::size_t sp, std::size_t ep) {
        sp_orig = sp;
        ep_orig = ep;
      }
    } state{*rank_};

    std::size_t run_sp = t_sp, run_ep = t_ep;
    PreprocessILCP<IlcpState>{state}(run_sp, run_ep);

    MarkedReported mr(n_doc_);
    const auto& select = *select_;

    auto get_doc = [this, &select, &state](std::size_t i) {
      const std::size_t head = static_cast<std::size_t>(select(i + 1));
      return get_doc_(std::max(state.sp_orig, head));
    };

    auto stop_pred = [this, t_m](std::size_t k, std::size_t /*bp*/) {
      return static_cast<std::size_t>((*run_values_)[k]) >= t_m;
    };

    auto report = [this, &mr, &t_report, &select, &state](std::size_t i, std::size_t doc) {
      auto report_dedup = [&mr, &t_report](auto d) {
        const auto dd = static_cast<std::size_t>(d);
        if (!mr(0, dd)) {
          mr.mark(dd);
          t_report(dd);
        }
      };
      // See IlcpLikeLeanCore::findDocs: CILCP★ merges by document too, so a visited
      // run must be fanned out even when its head doc is already reported.
      report_dedup(doc);

      const std::size_t head = static_cast<std::size_t>(select(i + 1));
      const std::size_t run_start = std::max(state.sp_orig, head);
      const std::size_t next_head = (i + 1 < n_runs_) ? static_cast<std::size_t>(select(i + 2)) : run_heads_->size();
      const std::size_t run_end = std::min(state.ep_orig - 1, next_head - 1);

      const std::size_t b = run_start + 1;
      const std::size_t e = run_end + 1;
      if constexpr (kVariantS == IlcpFullVariant::CILCP) {
        // CILCP★ runs can be either single-doc (merged) or multi-doc (a
        // non-merged, hence ilcp-constant, ILCP run). The peek is sound for
        // both, by the argument spelled out in IlcpLikeLeanCore::findDocs.
        if (b < e && get_doc_(b) != doc)
          get_doc_(b, e, report_dedup);
      } else {
        if (b < e)
          get_doc_(b, e, report_dedup);
      }
    };

    constexpr bool kAlwaysReport = (kVariantS == IlcpFullVariant::CILCP);
    ListDocsRMQSchemeDepth<kAlwaysReport>(run_sp, run_ep, *rmq_, get_doc, stop_pred, mr, report);
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    written += rmq_ ? sdsl::serialize(*rmq_, out, child, "rmq") : sdsl::serialize_empty_object<TRMQ>(out, child, "rmq");
    written += run_heads_ ? sdsl::serialize(*run_heads_, out, child, "run_heads")
                          : sdsl::serialize_empty_object<TBvRunHeads>(out, child, "run_heads");
    written += rank_ ? sdsl::serialize(*rank_, out, child, "run_heads_rank")
                     : sdsl::serialize_empty_object<typename TBvRunHeads::rank_1_type>(out, child, "run_heads_rank");
    written += select_
                   ? sdsl::serialize(*select_, out, child, "run_heads_select")
                   : sdsl::serialize_empty_object<typename TBvRunHeads::select_1_type>(out, child, "run_heads_select");
    written += run_values_
                   ? sdsl::serialize(*run_values_, out, child, "run_values")
                   : sdsl::serialize_empty_object<TRunValues>(out, child, "run_values");
    sdsl::int_vector<64> n_doc_vec(1, n_doc_);
    written += sdsl::serialize(n_doc_vec, out, child, "n_doc");
    written += get_doc_.serialize(out, child, "get_doc");
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rmq_)
      append(r, "rmq", sdsl::size_in_bytes(*rmq_));
    if (run_heads_)
      append(r, "run_heads", sdsl::size_in_bytes(*run_heads_));
    if (run_values_)
      append(r, "run_values", sdsl::size_in_bytes(*run_values_));
    auto get_doc_report = get_doc_.GetSizeReport();
    for (const auto& [k, v] : get_doc_report)
      append(r, k, v);
    return r;
  }

  TGetDoc& get_doc_policy() {
    return get_doc_;
  }

 protected:
  using typename Base::TSource;

  // One namespace per family holds everything the family needs: the partition
  // (run_heads, rmq), shared with the -L core, and the run values, which only
  // this core reads. ILCP and ILCP-L share the kILCP runs; CILCP and CILCP-L
  // share the kCILCP ones, both built by internal::BuildCilcpRuns.
  static std::string_view TopKey() {
    return kVariantS == IlcpFullVariant::ILCP ? conf::kILCP : conf::kCILCP;
  }

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    using namespace dret::conf;
    std::string top(TopKey());
    key_rmq_ = t_keys[top][kRmq].get<std::string>();
    key_run_heads_ = t_keys[top][kRunHeads].get<std::string>();
    key_run_values_ = t_keys[top][kRunValues].get<std::string>();

    rmq_ = this->template loadItemPtr<TRMQ>(key_rmq_, t_source, true);
    run_heads_ = this->template loadItemPtr<TBvRunHeads>(key_run_heads_, t_source, true);
    rank_ = &this->template loadBVRank<TBvRunHeads>(key_run_heads_, t_source, true).get();
    select_ = &this->template loadBVSelect<TBvRunHeads>(key_run_heads_, t_source, true).get();
    run_values_ = this->template loadItemPtr<TRunValues>(key_run_values_, t_source, true);

    n_runs_ = (*rank_)(run_heads_->size());

    const auto key_n_doc = t_keys[kRmqNDoc].get<std::string>();
    const auto* n_doc_vec = this->template loadItemPtr<sdsl::int_vector<64>>(key_n_doc, t_source, true);
    n_doc_ = static_cast<std::size_t>((*n_doc_vec)[0]);
  }

  std::string key_rmq_;
  std::string key_run_heads_;
  std::string key_run_values_;

  const TRMQ* rmq_ = nullptr;
  const TBvRunHeads* run_heads_ = nullptr;
  const typename TBvRunHeads::rank_1_type* rank_ = nullptr;
  const typename TBvRunHeads::select_1_type* select_ = nullptr;
  const TRunValues* run_values_ = nullptr;
  std::size_t n_doc_ = 0;
  std::size_t n_runs_ = 0;
  TGetDoc get_doc_;
};

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TBvRunHeads = sdsl::sd_vector<>,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>,
          typename TGetDoc = GetDocDA<TStorage, t_width>,
          typename TRunValues = sdsl::dac_vector<>>
using IlcpCore = IlcpLikeFullCore<IlcpFullVariant::ILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc, TRunValues>;

template <typename TStorage = GenericStorage,
          uint8_t t_width = 8,
          typename TBvRunHeads = sdsl::sd_vector<>,
          typename TRMQ = sdsl::rmq_succinct_sct<true>,
          typename TBvDocEnds = sdsl::sd_vector<>,
          typename TGetDoc = GetDocDA<TStorage, t_width>,
          typename TRunValues = sdsl::dac_vector<>>
using CilcpCore = IlcpLikeFullCore<IlcpFullVariant::CILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc, TRunValues>;

//~~~~~~~
// DocListIdxRMQ


template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TCountIdx = sri::RIndexCount<TStorage, TAlphabet>,
          typename TCore = SadaLCore<TStorage, TAlphabet::int_width>>
class DocListIdxRMQ : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using Base = DocListIndexExtStorage<TStorage, TAlphabet>;
  using Core = TCore;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;

  explicit DocListIdxRMQ(const TStorage& t_storage) : Base(t_storage), count_idx_(t_storage), core_(t_storage) {}

  DocListIdxRMQ(const TStorage& t_storage, const TCountIdx& t_count_idx)
      : Base(t_storage), count_idx_(t_count_idx), core_(t_storage) {}

  DocListIdxRMQ(const TStorage& t_storage, const TCore& t_core)
      : Base(t_storage), count_idx_(t_storage), core_(t_core) {}

  DocListIdxRMQ(const TStorage& t_storage, const TCountIdx& t_count_idx, const TCore& t_core)
      : Base(t_storage), count_idx_(t_count_idx), core_(t_core) {}

  DocListIdxRMQ() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    auto [sp, ep] = count_idx_.Count(t_pattern);
    if (sp >= ep)
      return;
    core_.findDocs(sp, ep, t_pattern.size(), t_report);
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
    for (const auto& [k, v] : core_report)
      append(r, k, v);
    return r;
  }

  const TCountIdx& count_idx() const {
    return count_idx_;
  }

  const TCore& core() const {
    return core_;
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

// Build the CILCP* run partition of Definition 1 in Cobas, Makinen and Rossi
// (SPIRE 2020) over the backward-ILCP array: identify the ILCP runs, flag each
// as single- or multi-document, then greedily merge consecutive single-document
// runs that share a document, keeping min(ILCP) over what is merged.
//
// Multi-document ILCP runs stand alone. Splitting one gains nothing -- the
// pieces keep the same value, since an ILCP run is value-uniform -- and costs a
// run head; merging one into a neighbour only lowers that neighbour's stored
// value and so weakens pruning. Either way the document-aware rule must leave
// them intact.
//
// Shared by CILCP and CILCP-S so the two cores partition IDENTICALLY by
// construction. They differ in exactly one thing: CILCP-S also stores the run
// values, which buys it the value-based stop; CILCP omits them and recurses
// unconditionally instead (see IlcpLikeLeanCore::findDocs).
inline void BuildCilcpRuns(const sdsl::int_vector<>& t_ilcp,
                           const sdsl::int_vector<>& t_da,
                           sdsl::bit_vector& t_run_heads,
                           std::vector<std::size_t>& t_run_values) {
  const std::size_t n = t_ilcp.size();
  t_run_heads = sdsl::bit_vector(n, 0);
  t_run_values.clear();
  if (n == 0) return;

  // Phase A -- the ILCP runs, each tagged with its document or kMultiDoc.
  constexpr std::size_t kMultiDoc = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> starts, values, docs;
  starts.push_back(0);
  values.push_back(t_ilcp[0]);
  std::size_t curr_doc = static_cast<std::size_t>(t_da[0]);
  bool curr_single = true;
  for (std::size_t i = 1; i < n; ++i) {
    if (t_ilcp[i] == t_ilcp[i - 1]) {
      if (static_cast<std::size_t>(t_da[i]) != curr_doc) curr_single = false;
    } else {
      docs.push_back(curr_single ? curr_doc : kMultiDoc);
      starts.push_back(i);
      values.push_back(t_ilcp[i]);
      curr_doc = static_cast<std::size_t>(t_da[i]);
      curr_single = true;
    }
  }
  docs.push_back(curr_single ? curr_doc : kMultiDoc);

  // Phase B -- greedy merge of consecutive single-document runs sharing a doc.
  const std::size_t n_runs = starts.size();
  std::size_t i = 0;
  while (i < n_runs) {
    t_run_heads[starts[i]] = 1;
    std::size_t v = values[i];
    if (docs[i] != kMultiDoc) {
      const std::size_t d = docs[i];
      std::size_t j = i + 1;
      while (j < n_runs && docs[j] == d) {
        if (values[j] < v) v = values[j];
        ++j;
      }
      t_run_values.push_back(v);
      i = j;
    } else {
      t_run_values.push_back(v);
      ++i;
    }
  }
}

}  // namespace internal

// Helper: pack a std::vector<std::size_t> of values into the requested
// TContainer and persist it under t_key. Used by ILCP-S / CILCP-S to
// store run_values (per-run min(VILCP)) and by SADA-S to store prev_doc
// (per SA-position previous occurrence). Compressed containers
// (sdsl::dac_vector, sdsl::vlc_vector, ...) compress at construction;
// sdsl::int_vector<> is bit-compressed explicitly to ceil(log2(max))
// bits per entry. Query-time access is operator[] in every case.
template <typename TContainer>
inline void StorePackedValues(Config& t_config,
                              const std::string& t_key,
                              const std::vector<std::size_t>& t_values) {
  if constexpr (std::is_same_v<TContainer, sdsl::int_vector<>>) {
    sdsl::int_vector<> packed(t_values.size());
    for (std::size_t i = 0; i < t_values.size(); ++i)
      packed[i] = t_values[i];
    sdsl::util::bit_compress(packed);
    sdsl::store_to_cache(packed, t_key, t_config, true);
  } else {
    std::vector<uint64_t> values(t_values.begin(), t_values.end());
    TContainer packed(values);
    sdsl::store_to_cache(packed, t_key, t_config, true);
  }
}

// Backwards-compatible alias for the run_values use site. Default
// TRunValues = sdsl::dac_vector<> matches the IlcpLikeFullCore default.
template <typename TRunValues = sdsl::dac_vector<>>
inline void StoreRunValues(Config& t_config,
                           const std::string& t_key,
                           const std::vector<std::size_t>& t_run_values) {
  StorePackedValues<TRunValues>(t_config, t_key, t_run_values);
}

// SADA construction: prev_doc + RMinQ.
template <typename TStorage, uint8_t t_width, typename TRMQ, typename TBvDocEnds, typename TGetDoc>
void construct(SadaLCore<TStorage, t_width, TRMQ, TBvDocEnds, TGetDoc>& t_core, Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kSADA][kRmq].get<std::string>();
  if (!sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)) {
    auto event = sdsl::memory_monitor::event(key_rmq);

    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

    std::size_t n_doc = internal::ReadNDoc(t_config);

    // prev_doc[i] = the previous SA position holding DA[i], or 0 when there is
    // none. Using 0 as "no previous occurrence" is only safe because SA position
    // 0 is the sentinel suffix, which no pattern range ever contains: every
    // query has sp >= 1, so 0 < sp acts as -infinity in the stop test
    // prev_doc[k] >= sp. If the sentinel handling ever changes so that position
    // 0 can fall inside a range, this needs a real out-of-band value.
    sdsl::int_vector<> prev_doc(da.size(), 0, sdsl::bits::hi(da.size()) + 1);
    std::vector<std::size_t> last_occ(n_doc + 2, 0);
    for (std::size_t i = 0; i < da.size(); ++i) {
      std::size_t doc = da[i];
      if (doc >= last_occ.size())
        last_occ.resize(doc + 1, 0);
      prev_doc[i] = last_occ[doc];
      last_occ[doc] = i;
    }

    TRMQ rmq(&prev_doc);
    sdsl::store_to_cache(rmq, key_rmq, t_config, true);
  }

  construct(t_core.get_doc_policy(), t_config);
}

// SADA-S construction: same prev_doc + RMinQ as SADA, plus a persisted
// prev_doc array (used by the depth-based recursion-stop at query time).
// TPrevDoc controls how that array is encoded on disk — int_vector is the
// natural default (random SA positions don't compress much under DAC/VLC).
template <typename TStorage, uint8_t t_width, typename TRMQ, typename TBvDocEnds, typename TGetDoc, typename TPrevDoc>
void construct(SadaCore<TStorage, t_width, TRMQ, TBvDocEnds, TGetDoc, TPrevDoc>& t_core, Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kSADA][kRmq].get<std::string>();
  auto key_prev_doc = t_config.keys[kSADA][kPrevDoc].get<std::string>();
  if (!sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      || !sdsl::cache_file_exists<TPrevDoc>(key_prev_doc, t_config)) {
    auto event = sdsl::memory_monitor::event(key_prev_doc);

    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

    std::size_t n_doc = internal::ReadNDoc(t_config);

    sdsl::int_vector<> prev_doc(da.size(), 0, sdsl::bits::hi(da.size()) + 1);
    std::vector<std::size_t> last_occ(n_doc + 2, 0);
    for (std::size_t i = 0; i < da.size(); ++i) {
      std::size_t doc = da[i];
      if (doc >= last_occ.size())
        last_occ.resize(doc + 1, 0);
      prev_doc[i] = last_occ[doc];
      last_occ[doc] = i;
    }

    if (!sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)) {
      TRMQ rmq(&prev_doc);
      sdsl::store_to_cache(rmq, key_rmq, t_config, true);
    }
    // Pack into the requested TPrevDoc container. For TPrevDoc =
    // sdsl::int_vector<> (the SadaCore default) StorePackedValues does
    // bit_compress; for DAC / VLC / similar it constructs from a
    // std::vector<uint64_t>.
    std::vector<std::size_t> prev_doc_values(prev_doc.size());
    for (std::size_t i = 0; i < prev_doc.size(); ++i)
      prev_doc_values[i] = prev_doc[i];
    StorePackedValues<TPrevDoc>(t_config, key_prev_doc, prev_doc_values);
  }

  construct(t_core.get_doc_policy(), t_config);
}

// ILCP construction: plain RLE on backward-ILCP.
template <typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds,
          typename TGetDoc>
void construct(IlcpLikeLeanCore<IlcpLeanVariant::ILCP_L, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc>& t_core,
               Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kILCP][kRmq].get<std::string>();
  auto key_run_heads = t_config.keys[kILCP][kRunHeads].get<std::string>();
  if (!sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      || !sdsl::cache_file_exists<TBvRunHeads>(key_run_heads, t_config)) {
    auto event = sdsl::memory_monitor::event(key_rmq);

    std::size_t n_doc = internal::ReadNDoc(t_config);
    auto ilcp = internal::LoadOrComputeIlcp<t_width>(t_config, n_doc);

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

  construct(t_core.get_doc_policy(), t_config);
}

// CILCP construction: the CMR20 CILCP* partition, run values not stored.
template <typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds,
          typename TGetDoc>
void construct(IlcpLikeLeanCore<IlcpLeanVariant::CILCP_L, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc>& t_core,
               Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kCILCP][kRmq].get<std::string>();
  auto key_run_heads = t_config.keys[kCILCP][kRunHeads].get<std::string>();
  if (!sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      || !sdsl::cache_file_exists<TBvRunHeads>(key_run_heads, t_config)) {
    auto event = sdsl::memory_monitor::event(key_rmq);

    std::size_t n_doc = internal::ReadNDoc(t_config);
    auto ilcp = internal::LoadOrComputeIlcp<t_width>(t_config, n_doc);

    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

    // Same partition as CILCP-S, by construction. CILCP builds the RMQ over the
    // run values but does NOT store them: without frequencies there is nothing
    // the value-based stop is needed for, and dropping the array is the whole
    // space saving of this core. The traversal pays for it by recursing
    // unconditionally (Proposition: the marker-based stop is unsound once runs
    // are merged by document).
    sdsl::bit_vector run_heads;
    std::vector<std::size_t> run_values;
    internal::BuildCilcpRuns(ilcp, da, run_heads, run_values);

    internal::StoreRunHeadsAndRMQ<TBvRunHeads, TRMQ>(
        t_config, key_run_heads, key_rmq, std::move(run_heads), run_values);
  }

  construct(t_core.get_doc_policy(), t_config);
}

// StorePackedValues / StoreRunValues are defined earlier in the file
// (right before construct(SadaLCore, ...) — see above).

// ILCP-S construction: reuse the existing ILCP RLE (run_heads + rmq) and
// additionally persist run_values for the depth-based stop. If the ILCP
// caches already exist (e.g., built by the existing ILCP construct on a
// prior run), we only recompute the values; otherwise we build both.
template <typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds,
          typename TGetDoc,
          typename TRunValues>
void construct(IlcpLikeFullCore<IlcpFullVariant::ILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc, TRunValues>& t_core,
               Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kILCP][kRmq].get<std::string>();
  auto key_run_heads = t_config.keys[kILCP][kRunHeads].get<std::string>();
  auto key_run_values = t_config.keys[kILCP][kRunValues].get<std::string>();

  const bool rle_missing =
      !sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      || !sdsl::cache_file_exists<TBvRunHeads>(key_run_heads, t_config);
  const bool values_missing = !sdsl::cache_file_exists<TRunValues>(key_run_values, t_config);

  if (rle_missing || values_missing) {
    auto event = sdsl::memory_monitor::event(key_run_values);

    std::size_t n_doc = internal::ReadNDoc(t_config);
    auto ilcp = internal::LoadOrComputeIlcp<t_width>(t_config, n_doc);

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

    if (rle_missing) {
      internal::StoreRunHeadsAndRMQ<TBvRunHeads, TRMQ>(
          t_config, key_run_heads, key_rmq, std::move(run_heads), run_values);
    }
    if (values_missing) {
      StoreRunValues<TRunValues>(t_config, key_run_values, run_values);
    }
  }

  construct(t_core.get_doc_policy(), t_config);
}

// CILCP-S construction: reuse CILCP's CMR20 partition (run_heads + rmq under
// kCILCP, built by internal::BuildCilcpRuns) and additionally persist run_values
// for the value-based stop. Mirrors ILCP-S over ILCP. If the CILCP caches already
// exist we only compute the values; otherwise we build both. Storing one copy of
// a partition the two cores share by construction removes the duplicate on disk,
// the second ILCP+DA scan, and the chance of the two copies drifting apart.
template <typename TStorage,
          uint8_t t_width,
          typename TBvRunHeads,
          typename TRMQ,
          typename TBvDocEnds,
          typename TGetDoc,
          typename TRunValues>
void construct(IlcpLikeFullCore<IlcpFullVariant::CILCP, TStorage, t_width, TBvRunHeads, TRMQ, TBvDocEnds, TGetDoc, TRunValues>& t_core,
               Config& t_config) {
  using namespace dret::conf;
  internal::EnsureBasicStructures<t_width, TBvDocEnds>(t_config);

  auto key_rmq = t_config.keys[kCILCP][kRmq].get<std::string>();
  auto key_run_heads = t_config.keys[kCILCP][kRunHeads].get<std::string>();
  auto key_run_values = t_config.keys[kCILCP][kRunValues].get<std::string>();

  const bool rle_missing =
      !sdsl::cache_file_exists<TRMQ>(key_rmq, t_config)
      || !sdsl::cache_file_exists<TBvRunHeads>(key_run_heads, t_config);
  const bool values_missing = !sdsl::cache_file_exists<TRunValues>(key_run_values, t_config);

  if (rle_missing || values_missing) {
    auto event = sdsl::memory_monitor::event(key_run_values);

    std::size_t n_doc = internal::ReadNDoc(t_config);
    auto ilcp = internal::LoadOrComputeIlcp<t_width>(t_config, n_doc);

    sdsl::int_vector<> da;
    sdsl::load_from_cache(da, t_config.keys[kDA].get<std::string>(), t_config, true);

    // Identical partition to CILCP, by construction -- the two cores differ only
    // in that CILCP-S also stores the run values below, which is what enables
    // its value-based stop.
    sdsl::bit_vector run_heads;
    std::vector<std::size_t> run_values;
    internal::BuildCilcpRuns(ilcp, da, run_heads, run_values);

    if (rle_missing) {
      internal::StoreRunHeadsAndRMQ<TBvRunHeads, TRMQ>(
          t_config, key_run_heads, key_rmq, std::move(run_heads), run_values);
    }
    if (values_missing) {
      StoreRunValues<TRunValues>(t_config, key_run_values, run_values);
    }
  }

  construct(t_core.get_doc_policy(), t_config);
}

// Wiring: construct the whole DocListIdxRMQ (count_idx + core).
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TCore>
void construct(DocListIdxRMQ<TStorage, TAlphabet, TCountIdx, TCore>& t_index, Config& t_config) {
  auto count_idx = t_index.count_idx();
  construct(count_idx, t_config.data_path, t_config);

  TCore core(t_index.core());
  construct(core, t_config);

  t_index.load(t_config);
}

}  // namespace rmq
}  // namespace dret
