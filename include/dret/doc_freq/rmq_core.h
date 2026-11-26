//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Document frequencies with the RMQ cores (SADA, ILCP, CILCP), as the 2020 frequency code computed
// them: the leftmost and the rightmost position of each document in the range, and a frequency oracle
// (oracle.h) that turns the two into the document's frequency.
//
// The leftmost positions come from the published listing core itself (doc_list/doc_list_rmq.h): its
// marking traversal reports each document from its first occurrence in the range, and hands the
// position to a sink that takes it (ReportPosition, sink.h). The rightmost positions are the leftmost
// ones of the reversed problem -- position i >= 1 of the suffix array is position n - i of the reversed
// one, and position 0, the terminator's, which no pattern range holds, stays first, as the listing
// constructions assume -- so they come from a second listing core of the same kind, built over the reversed document
// array (and, for the ILCP family, the reversed forward ILCP array: the backward ILCP of the reversed
// order) by the listing core's own construction, and reading the documents through a reversed view of
// the first core's backend. Its structures are cached under the listing keys with a "_rev" suffix
// (RightKeys).
//
// A core of this kind is TCoreOf<TGetDoc>: the listing core with its backend as the template argument.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "dret/collection.h"
#include "dret/config.h"
#include "dret/doc_freq/frequencies.h"
#include "dret/doc_freq/oracle.h"
#include "dret/doc_list/distinct_docs.h"
#include "dret/doc_list/doc_list_rmq.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret::freq {

namespace internal {

// The leftmost position of each document a listing core reports (ReportPosition).
class ExtremesSink {
 public:
  explicit ExtremesSink(std::size_t t_n_doc) : seen_(t_n_doc) {}

  void PositionAt(std::size_t t_doc, std::size_t t_pos) {
    if (seen_.Insert(t_doc)) found_.emplace_back(t_doc, t_pos);
  }

  // Not called: the RMQ cores report through ReportPosition. Here for the sink concept.
  void Position(std::size_t t_doc) { PositionAt(t_doc, 0); }

  void Listed(std::size_t t_doc) { Position(t_doc); }

  bool Contains(std::size_t t_doc) const { return seen_.Contains(t_doc); }

  static constexpr bool Done() { return false; }

  // (document, position), sorted by document.
  std::vector<std::pair<std::size_t, std::size_t>>& Sorted() {
    std::sort(found_.begin(), found_.end());
    return found_;
  }

 private:
  DistinctDocs seen_;
  std::vector<std::pair<std::size_t, std::size_t>> found_;
};

}  // namespace internal

// The document of position j >= 1 of the reversed suffix array: the backend's at n - j (position 0 is
// the terminator's in both). A view of the backend the left core holds, bound before each query; it
// stores and builds nothing.
template <typename TBackend, typename TStorage = GenericStorage, uint8_t t_width = 8>
class ReversedView : public FieldsComponent<ReversedView<TBackend, TStorage, t_width>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<ReversedView, TStorage, t_width>;

  explicit ReversedView(const TStorage& t_storage) : Base(t_storage) {}

  ReversedView() = default;

  void Bind(const TBackend* t_backend, std::size_t t_n) const {
    backend_ = t_backend;
    n_ = t_n;
  }

  std::size_t operator()(std::size_t t_j) const { return (*backend_)(Original(t_j)); }

  template <typename TReport>
  void operator()(std::size_t t_b, std::size_t t_e, TReport& t_report) const {
    for (std::size_t j = t_b; j < t_e; ++j) t_report((*this)(j));
  }

  // In order of the reversed positions: the backend's range read forward, then walked backwards.
  template <typename F>
  void ExpandUntil(std::size_t t_b, std::size_t t_e, F&& t_f) const {
    if (t_b >= t_e) return;
    if (t_b == 0) {  // the terminator's position, the same in both
      if (!t_f((*backend_)(0))) return;
      ++t_b;
      if (t_b >= t_e) return;
    }
    std::vector<std::size_t> docs;
    docs.reserve(t_e - t_b);
    backend_->ExpandUntil(Original(t_e - 1), Original(t_b) + 1, [&docs](std::size_t d) {
      docs.push_back(d);
      return true;
    });
    for (auto it = docs.rbegin(); it != docs.rend(); ++it)
      if (!t_f(*it)) return;
  }

  std::size_t size() const { return n_; }

  // The position of the suffix array that position t_j of the reversed one stands for.
  std::size_t Original(std::size_t t_j) const { return t_j == 0 ? 0 : n_ - t_j; }

 protected:
  friend Base;

  template <typename TSelf>
  static auto Fields(TSelf&) {
    return std::tuple{};
  }

  mutable const TBackend* backend_ = nullptr;
  mutable std::size_t n_ = 0;
};

template <typename TBackend, typename TStorage, uint8_t t_width>
void construct(ReversedView<TBackend, TStorage, t_width>&, Config&) {}

// The keys of the reversed structures: the document array, the ILCP array and every SADA / ILCP /
// CILCP structure, with a "_rev" suffix.
inline JSON RightKeys(const JSON& t_keys) {
  JSON keys = t_keys;
  for (auto name : {conf::kDA, conf::kIlcpArray}) keys[name] = keys[name].get<std::string>() + "_rev";
  for (auto core : {conf::kSADA, conf::kILCP, conf::kCILCP})
    for (auto& [k, v] : keys[core].items()) v = v.get<std::string>() + "_rev";
  return keys;
}

template <typename T>
struct IsIlcpFamily : std::false_type {};
template <rmq::IlcpPartition kVariant, typename TStorage, uint8_t t_width, typename... Ts>
struct IsIlcpFamily<rmq::IlcpLikeCore<kVariant, TStorage, t_width, Ts...>> : std::true_type {};
// The reversed forward ILCP array is built only for the ILCP family: a trait that missed it would leave
// the reversed core to compute a backward ILCP over the reversed document array, which is meaningless.
static_assert(IsIlcpFamily<rmq::IlcpCore<>>::value && IsIlcpFamily<rmq::CilcpCore<>>::value &&
              !IsIlcpFamily<rmq::SadaCore<>>::value);

template <typename TStorage, uint8_t t_width, typename TBackend, template <typename> class TCoreOf, typename TOracle>
class RmqFreqCore : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;
  using Left = TCoreOf<TBackend>;
  using Right = TCoreOf<ReversedView<TBackend, TStorage, t_width>>;

  // t_block_size: the sample spacing of the grammars the backend and the oracle keep, when they keep one
  // (backend::DiffSA, DocIsaOracle); the others ignore it.
  explicit RmqFreqCore(const TStorage& t_storage, uint32_t t_block_size = 256)
      : Base(t_storage), left_(MakeLeft(t_storage, t_block_size)), right_(t_storage),
        oracle_(MakeOracle(t_storage, t_block_size)) {}

  RmqFreqCore() = default;

  // A listing sink gets the listing core's answer; a frequency sink, the documents with their frequencies.
  template <typename TSink>
  void Query(std::size_t t_sp, std::size_t t_ep, std::size_t t_m, TSink& t_sink) const {
    if constexpr (!requires { t_sink.Weighted(std::size_t{}, std::size_t{}); }) {
      left_.Query(t_sp, t_ep, t_m, t_sink);
    } else {
      QueryFreq(t_sp, t_ep, t_m, t_sink);
    }
  }

  template <typename TSink>
  void QueryFreq(std::size_t t_sp, std::size_t t_ep, std::size_t t_m, TSink& t_sink) const {
    if (t_sp >= t_ep) return;
    const auto& backend = left_.get_doc_policy();
    const std::size_t n = backend.size();

    internal::ExtremesSink leftmost(left_.n_doc());
    left_.Query(t_sp, t_ep, t_m, leftmost);
    right_.get_doc_policy().Bind(&backend, n);
    internal::ExtremesSink rightmost(left_.n_doc());
    // [sp, ep) with sp >= 1 is [n - ep + 1, n - sp + 1) reversed.
    right_.Query(n - t_ep + 1, n - t_sp + 1, t_m, rightmost);

    const auto& l = leftmost.Sorted();
    const auto& r = rightmost.Sorted();
    for (std::size_t i = 0; i < l.size() && i < r.size(); ++i) {
      const auto d = l[i].first;
      t_sink.Weighted(d, oracle_.Freq(backend, d, l[i].second, n - r[i].second));
    }
  }

  std::size_t n_doc() const { return left_.n_doc(); }

  void load(Config t_config) override {
    left_.load(t_config);
    auto right_config = t_config;
    right_config.keys = RightKeys(t_config.keys);
    right_.load(right_config);
    oracle_.load(t_config);
  }

  void load(std::istream& t_in, const JSON& t_keys) override {
    left_.load(t_in, t_keys);
    right_.load(t_in, RightKeys(t_keys));
    oracle_.load(t_in, t_keys);
  }

  using Base::load;

  size_type serialize(std::ostream& t_out, sdsl::structure_tree_node* t_v, const std::string& t_name) const override {
    auto child = sdsl::structure_tree::add_child(t_v, t_name, sdsl::util::class_name(*this));
    size_type written = left_.serialize(t_out, child, "left");
    written += right_.serialize(t_out, child, "right");
    written += oracle_.serialize(t_out, child, "oracle");
    sdsl::structure_tree::add_size(child, written);
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    for (const auto& [name, bytes] : left_.GetSizeReport()) append(r, name, bytes);
    for (const auto& [name, bytes] : right_.GetSizeReport()) append(r, "rev_" + name, bytes);
    for (const auto& [name, bytes] : oracle_.GetSizeReport()) append(r, "oracle_" + name, bytes);
    return r;
  }

  static Left MakeLeft(const TStorage& t_storage, uint32_t t_block_size) {
    if constexpr (std::is_constructible_v<TBackend, const TStorage&, uint32_t, float>)
      return Left(t_storage, t_block_size, 0.0f);
    else
      return Left(t_storage);
  }

  static TOracle MakeOracle(const TStorage& t_storage, uint32_t t_block_size) {
    if constexpr (std::is_constructible_v<TOracle, const TStorage&, uint32_t>)
      return TOracle(t_storage, t_block_size);
    else
      return TOracle(t_storage);
  }

  const Left& left() const { return left_; }
  const Right& right() const { return right_; }
  const TOracle& oracle() const { return oracle_; }

 protected:
  Left left_;
  Right right_;
  TOracle oracle_;
};

// The listing core's structures; the reversed document array (and the reversed forward ILCP array for
// the ILCP family) under the "_rev" keys; the reversed core's structures, by the listing construction
// over them; the oracle's.
template <typename TStorage, uint8_t t_width, typename TBackend, template <typename> class TCoreOf, typename TOracle>
void construct(RmqFreqCore<TStorage, t_width, TBackend, TCoreOf, TOracle>& t_core, Config& t_config) {
  using Core = RmqFreqCore<TStorage, t_width, TBackend, TCoreOf, TOracle>;
  {
    auto left = t_core.left();
    construct(left, t_config);
  }

  auto right_config = t_config;
  right_config.keys = RightKeys(t_config.keys);
  const auto key_da = t_config.keys[conf::kDA].get<std::string>();
  const auto key_da_rev = right_config.keys[conf::kDA].get<std::string>();
  sdsl::int_vector<> da;
  if (!sdsl::cache_file_exists<sdsl::int_vector<>>(key_da_rev, t_config)) {
    sdsl::load_from_cache(da, key_da, t_config, true);
    // Position 0, the terminator's, stays first; position i >= 1 goes to n - i.
    sdsl::int_vector<> rev(da.size(), 0, da.width());
    rev[0] = da[0];
    for (std::size_t i = 1; i < da.size(); ++i) rev[da.size() - i] = da[i];
    sdsl::store_to_cache(rev, key_da_rev, t_config, true);
  }
  if constexpr (IsIlcpFamily<typename Core::Left>::value) {
    const auto key_ilcp_rev = right_config.keys[conf::kIlcpArray].get<std::string>();
    if (!sdsl::cache_file_exists<sdsl::int_vector<>>(key_ilcp_rev, t_config)) {
      if (da.empty()) sdsl::load_from_cache(da, key_da, t_config, true);
      const auto bwd = rmq::internal::LoadOrComputeIlcp<t_width>(t_config, ReadNDoc<t_width>(t_config));
      // The forward ILCP at i is the backward ILCP at the next position of the same document (0 when
      // there is none); stored reversed.
      sdsl::int_vector<> rev(bwd.size(), 0, bwd.width());
      constexpr std::size_t kNone = static_cast<std::size_t>(-1);
      std::vector<std::size_t> next;
      for (std::size_t i = da.size(); i-- > 0;) {
        const std::size_t d = da[i];
        if (d >= next.size()) next.resize(d + 1, kNone);
        if (i > 0) rev[da.size() - i] = next[d] != kNone ? bwd[next[d]] : 0;
        next[d] = i;
      }
      sdsl::store_to_cache(rev, key_ilcp_rev, t_config, true);
    }
  }
  {
    auto right = t_core.right();
    construct(right, right_config);
  }
  {
    auto oracle = t_core.oracle();
    construct(oracle, t_config);
  }
}

// Every document of the range with its frequency (Weighted).
template <typename TStorage, uint8_t t_width, typename TBackend, template <typename> class TCoreOf, typename TOracle>
struct ComputesFrequencies<RmqFreqCore<TStorage, t_width, TBackend, TCoreOf, TOracle>> : std::true_type {};

// The published listing cores with their backend as the only parameter, for RmqFreqCore.
template <typename TStorage = GenericStorage, uint8_t t_width = 8>
struct Cores {
  template <typename TGetDoc>
  using Sada = rmq::SadaCore<TStorage, t_width, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
  template <typename TGetDoc>
  using Ilcp = rmq::IlcpCore<TStorage, t_width, sdsl::sd_vector<>, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
  template <typename TGetDoc>
  using Cilcp = rmq::CilcpCore<TStorage, t_width, sdsl::sd_vector<>, sdsl::rmq_succinct_sct<true>, sdsl::sd_vector<>, TGetDoc>;
};

// An RMQ frequency index: the r-index for the range, then RmqFreqCore.
template <typename TStorage,
          typename TAlphabet,
          typename TCountIdx,
          typename TBackend,
          template <typename> class TCoreOf,
          typename TOracle>
using DocFreqIdxRMQ =
    DocIndex<TStorage, TAlphabet, TCountIdx, RmqFreqCore<TStorage, TAlphabet::int_width, TBackend, TCoreOf, TOracle>>;

}  // namespace dret::freq
