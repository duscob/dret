//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// The listing index: a range finder and a core. Every listing index is a DocListIdx:
//
//   range finder  [sp, ep) of the pattern: the r-index count (TCountIdx), or the core itself when it
//                 finds the range (a backend that counts, the RLCSA: backend::CountsPattern)
//   core          the documents of [sp, ep), into a sink (dret/sink.h):
//                   BruteCore        a backend expands the whole range           (brute_core.h)
//                   SampledTreeCore  the stored lists of the sampled nodes inside
//                                    the range, a backend for the rest (GCDA, PDL) (sampled_tree_core.h)
//                   the RMQ cores    a traversal of an RMQ partition             (doc_list_rmq.h)
//                 or the documents of the pattern's occurrences, found by the core itself (LocateCore,
//                 Brute-RI / Brute-SRI over the r-index / sr-index locate), which then needs no range.
//
// A core is an IndexBaseWithExternalStorage (load, serialize) with
//
//   Query(sp, ep, m, sink)   the documents of [sp, ep); m is the pattern length    (a range core)
//   Query(pattern, sink)     the documents of the pattern's occurrences            (a pattern core)
//   Count(pattern)           optional, a range core: [sp, ep) in place of the r-index
//   n_doc()                  the number of documents, or 0 when not loaded (the listing bitmap then grows)
//   GetSizeReport()          its fields
//   kFrequencies             optional, true when its Query serves document frequencies (FreqSink): it
//                            hands the sink every position of the range once, or the counts of the parts
//                            it does not expand
//
// and a construct(core, config, ...) that builds its structures.
//
// The index holds what every family shares: the r-index (loaded, serialized, sized and built only when
// the core needs it), the listing Search, and the serialization order -- the r-index, then the core --,
// which is the layout every index had before this wrapper, so the files and the sizes are unchanged.
//

#pragma once

#include <concepts>
#include <cstddef>
#include <functional>
#include <istream>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>

#include <sdsl/io.hpp>

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"
#include "dret/index_base.h"
#include "dret/sink.h"
#include "dret/size_report.h"

#ifdef DRET_DOC_LIST_PROFILE
#include <chrono>

#include "dret/doc_list/search_profile.h"
#endif

namespace dret {

namespace internal {
// A sink, for the core checks.
struct AnySink {
  void Position(std::size_t) {}
  void Listed(std::size_t) {}
  static constexpr bool Done() { return false; }
};
}  // namespace internal

// The core finds the occurrences of the pattern itself.
template <typename TCore, typename TPattern>
concept PatternCore = requires(const TCore& t_core, const TPattern& t_pattern, internal::AnySink& t_sink) {
  t_core.Query(t_pattern, t_sink);
};

// The core finds the suffix-array range of the pattern itself.
template <typename TCore, typename TPattern>
concept CoreCountsPattern = requires(const TCore& t_core, const TPattern& t_pattern) {
  { t_core.Count(t_pattern) } -> std::convertible_to<std::pair<std::size_t, std::size_t>>;
};

// The core's Query serves document frequencies (dret/sink.h, FreqSink).
template <typename TCore>
concept FrequencyCore = requires {
  requires TCore::kFrequencies;
};

// The range finder of an index whose core needs none: constructible from anything, holds nothing.
struct NoCountIdx {
  NoCountIdx() = default;
  template <typename... TArgs>
  explicit NoCountIdx(TArgs&&...) {}
};

template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TCore>
class DocListIdx : public DocListIndexExtStorage<TStorage, TAlphabet> {
 public:
  using Base = DocListIndexExtStorage<TStorage, TAlphabet>;
  using typename Base::size_type;
  using typename Base::TDocId;
  using typename Base::TPattern;
  using Core = TCore;

  static constexpr bool kPatternCore = PatternCore<TCore, TPattern>;
  static constexpr bool kCoreCounts = CoreCountsPattern<TCore, TPattern>;
  // The r-index is part of the index only when the core needs a range and does not find it itself.
  static constexpr bool kUsesCountIdx = !kPatternCore && !kCoreCounts;

  using CountIdx = std::conditional_t<kUsesCountIdx, TCountIdx, NoCountIdx>;

  // The core's parameters follow the storage (e.g. block size and storing factor).
  template <typename... TArgs>
    requires std::constructible_from<TCore, const TStorage&, TArgs...>
  explicit DocListIdx(const TStorage& t_storage, TArgs&&... t_core_args)
      : Base(t_storage), count_idx_(t_storage), core_(t_storage, std::forward<TArgs>(t_core_args)...) {}

  // A core set up by the caller (e.g. the RMQ cores over a backend with its own parameters).
  DocListIdx(const TStorage& t_storage, const TCore& t_core) : Base(t_storage), count_idx_(t_storage), core_(t_core) {}

  DocListIdx() = default;

  void Search(const TPattern& t_pattern, const std::function<void(TDocId)>& t_report) const override {
    ListSink<std::function<void(TDocId)>> sink(core_.n_doc(), t_report);
    Query(t_pattern, sink);
  }

  // The query into any sink (dret/sink.h).
  template <typename TSink>
  void Query(const TPattern& t_pattern, TSink& t_sink) const {
    if constexpr (kPatternCore) {
      core_.Query(t_pattern, t_sink);
    } else {
#ifdef DRET_DOC_LIST_PROFILE
      const auto prof_t0 = std::chrono::steady_clock::now();
#endif
      const auto [sp, ep] = Count(t_pattern);
#ifdef DRET_DOC_LIST_PROFILE
      if constexpr (requires { TCore::kProfiled; }) {
        search_profile().ns_count +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - prof_t0).count();
      }
#endif
      if (sp >= ep) return;
      core_.Query(sp, ep, t_pattern.size(), t_sink);
    }
  }

  // The documents of the pattern with their frequencies: t_report(d, f) once per document, in no given
  // order. Only over a core that serves frequencies.
  template <typename TReport>
    requires FrequencyCore<TCore>
  void Frequencies(const TPattern& t_pattern, const TReport& t_report) const {
    FreqSink<TReport> sink(core_.n_doc(), t_report);
    Query(t_pattern, sink);
    sink.Finish();
  }

  // The suffix-array range of the pattern (a range core only).
  std::pair<std::size_t, std::size_t> Count(const TPattern& t_pattern) const
    requires(!kPatternCore)
  {
    if constexpr (kCoreCounts) {
      return core_.Count(t_pattern);
    } else {
      const auto [sp, ep] = count_idx_.Count(t_pattern);
      return {sp, ep};
    }
  }

  void load(Config t_config) override {
    if constexpr (kUsesCountIdx) count_idx_.load(t_config);
    core_.load(t_config);
  }

  void load(std::istream& in, const JSON& t_keys) override {
    if constexpr (kUsesCountIdx) count_idx_.load(in);
    core_.load(in, t_keys);
  }

  using Base::load;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    if constexpr (kUsesCountIdx) written += count_idx_.serialize(out, child, "count_idx");
    written += core_.serialize(out, child, "core");
    return written;
  }

  SizeReport GetSizeReport() const override {
    SizeReport r;
    if constexpr (kUsesCountIdx) append(r, "count_idx", sdsl::size_in_bytes(count_idx_));
    for (const auto& [name, bytes] : core_.GetSizeReport()) append(r, name, bytes);
    return r;
  }

  const CountIdx& count_idx() const { return count_idx_; }

  const TCore& core() const { return core_; }

  TCore& core() { return core_; }

 protected:
  CountIdx count_idx_;
  TCore core_;
};

// The core's structures (on a copy of the core, which keeps the index's parameters), then the r-index
// when the index uses it, then the load. Extra arguments go to the core's construct.
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TCore, typename... TArgs>
void construct(DocListIdx<TStorage, TAlphabet, TCountIdx, TCore>& t_index, Config& t_config, TArgs&&... t_args) {
  {
    auto core = t_index.core();
    construct(core, t_config, std::forward<TArgs>(t_args)...);
  }
  if constexpr (DocListIdx<TStorage, TAlphabet, TCountIdx, TCore>::kUsesCountIdx) {
    auto count_idx = t_index.count_idx();
    construct(count_idx, t_config.data_path, t_config);
  }
  t_index.load(t_config);
}

}  // namespace dret
