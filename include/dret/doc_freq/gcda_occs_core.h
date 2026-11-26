//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Document frequencies over the DA grammar with occurrence bitvectors (GCDA-ISAs and GCDA-WT-DA in the
// 2020 frequency code). Every rule stores two bitvectors over the documents (common/tf.h): which
// documents' first occurrence in the rule's expansion is in its right child, and which documents'
// last occurrence is. The range [sp, ep) is covered by rules (the top-level symbols inside it and
// the partial covers of the two at its ends); joined in a virtual binary tree, they give each
// document's leftmost and rightmost positions in the range, which a frequency oracle (oracle.h) turns
// into its frequency with the backend's help.
//
// Stored: the grammar's rules bit-packed, its top-level sequence, where each top-level symbol starts
// in the DA, and the two bitvectors per rule; then the backend and the oracle.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/sd_vector.hpp>

#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "dret/collection.h"
#include "dret/common/tf.h"
#include "dret/config.h"
#include "dret/doc_freq/frequencies.h"
#include "dret/doc_freq/oracle.h"
#include "dret/doc_list/doc_list_index.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"
#include "dret/slp/da_grammar.h"
#include "dret/slp/slp_components.h"

namespace dret::freq {

namespace internal {

// The rules that cover [sp, ...) of t_var's expansion, left to right.
template <typename TSLP, typename TReport>
void CoverFrom(const TSLP& t_slp, std::size_t t_var, std::size_t t_sp, const TReport& t_report) {
  if (t_sp == 0) return t_report(t_var);
  const auto children = t_slp[t_var];
  const auto left = t_slp.SpanLength(children.first);
  if (left <= t_sp) return CoverFrom(t_slp, children.second, t_sp - left, t_report);
  CoverFrom(t_slp, children.first, t_sp, t_report);
  t_report(children.second);
}

// The rules that cover [0, ep] of t_var's expansion, left to right.
template <typename TSLP, typename TReport>
void CoverTo(const TSLP& t_slp, std::size_t t_var, std::size_t t_ep, const TReport& t_report) {
  if (t_slp.SpanLength(t_var) <= t_ep + 1) return t_report(t_var);
  const auto children = t_slp[t_var];
  const auto left = t_slp.SpanLength(children.first);
  if (t_ep < left) return CoverTo(t_slp, children.first, t_ep, t_report);
  t_report(children.first);
  CoverTo(t_slp, children.second, t_ep - left, t_report);
}

// The rules that cover [sp, ep] of t_var's expansion, left to right.
template <typename TSLP, typename TReport>
void CoverIn(const TSLP& t_slp, std::size_t t_var, std::size_t t_sp, std::size_t t_ep, const TReport& t_report) {
  if (t_sp > t_ep) return;
  if (t_sp == 0 && t_slp.SpanLength(t_var) <= t_ep + 1) return t_report(t_var);
  const auto children = t_slp[t_var];
  const auto left = t_slp.SpanLength(children.first);
  if (t_ep < left) return CoverIn(t_slp, children.first, t_sp, t_ep, t_report);
  if (left <= t_sp) return CoverIn(t_slp, children.second, t_sp - left, t_ep - left, t_report);
  CoverFrom(t_slp, children.first, t_sp, t_report);
  CoverTo(t_slp, children.second, t_ep - left, t_report);
}

}  // namespace internal

template <typename TStorage, uint8_t t_width, typename TBackend, typename TOracle>
class GcdaOccsCore : public FieldsComponent<GcdaOccsCore<TStorage, t_width, TBackend, TOracle>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<GcdaOccsCore, TStorage, t_width>;
  using typename Base::size_type;
  using SLP = grammar::SLP<sdsl::int_vector<>, sdsl::int_vector<>>;
  using OccsBVs = std::vector<sdsl::bit_vector>;
  using RootHeads = sdsl::sd_vector<>;

  // t_block_size: the sample spacing of the backend's and the oracle's grammars, when they keep one.
  explicit GcdaOccsCore(const TStorage& t_storage, uint32_t t_block_size = 256)
      : Base(t_storage), backend_(MakePart<TBackend>(t_storage, t_block_size)), oracle_(MakePart<TOracle>(t_storage, t_block_size)) {}

  GcdaOccsCore() = default;

  // A frequency sink gets every document with its frequency; a listing one, every document.
  template <typename TSink>
  void Query(std::size_t t_sp, std::size_t t_ep, std::size_t /*t_m*/, TSink& t_sink) const {
    if (t_sp >= t_ep) return;
    const std::size_t ep = t_ep - 1;  // inclusive from here, as the covers are

    std::vector<std::size_t> cover;
    auto add = [&cover](std::size_t t_var) { cover.push_back(t_var); };
    auto l = (*rank_)(t_sp + 1) - 1;
    const auto r = (*rank_)(ep + 1) - 1;
    const std::size_t l_start = (*select_)(l + 1);
    if (l < r) {
      if (l_start < t_sp) {
        internal::CoverFrom(*slp_, (*roots_)[l], t_sp - l_start, add);
        ++l;
      }
      for (auto i = l; i < r; ++i) cover.push_back((*roots_)[i]);
      internal::CoverTo(*slp_, (*roots_)[r], ep - (*select_)(r + 1), add);
    } else {
      internal::CoverIn(*slp_, (*roots_)[l], t_sp - l_start, ep - l_start, add);
    }

    // The first and the last occurrence of each document, 1-based from sp; the terminator's id, past
    // the documents, never occurs in a pattern range.
    std::vector<std::pair<std::size_t, std::size_t>> first, last;
    auto report_first = [&first](std::size_t d, std::size_t pos) { first.emplace_back(d, pos); };
    auto report_last = [&last](std::size_t d, std::size_t pos) { last.emplace_back(d, pos); };
    FindAllFirstLastOccs<grammar::SLP<>>(*slp_, cover, *first_occs_, *last_occs_, report_first, report_last);

    std::sort(first.begin(), first.end());
    std::sort(last.begin(), last.end());
    for (std::size_t i = 0; i < first.size() && i < last.size(); ++i) {
      const auto d = first[i].first;
      if (d >= n_doc_) continue;
      if constexpr (requires { t_sink.Weighted(std::size_t{}, std::size_t{}); }) {
        t_sink.Weighted(d, oracle_.Freq(backend_, d, t_sp + first[i].second - 1, t_sp + last[i].second - 1));
      } else {
        t_sink.Position(d);
      }
    }
  }

  std::size_t n_doc() const { return n_doc_; }

  const TBackend& backend() const { return backend_; }
  const TOracle& oracle() const { return oracle_; }

 protected:
  friend Base;

  template <typename TPart>
  static TPart MakePart(const TStorage& t_storage, uint32_t t_block_size) {
    if constexpr (std::is_constructible_v<TPart, const TStorage&, uint32_t>)
      return TPart(t_storage, t_block_size);
    else
      return TPart(t_storage);
  }

  // The structures, in the order they are serialized and read back (dret/fields.h).
  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    using namespace conf;
    auto key = [](std::string_view t_name) { return [t_name](const JSON& k) { return KeyName(k, t_name); }; };
    return std::tuple{Item("slp", t_self.slp_, key(kDaGrammarPacked)),
                      Item("roots", t_self.roots_, key(kDaGrammarSeq)),
                      Item("root_heads", t_self.root_heads_, key(kDaGrammarRootHeads)),
                      Rank<RootHeads>("root_heads_rank", t_self.rank_, key(kDaGrammarRootHeads)),
                      Select<RootHeads>("root_heads_select", t_self.select_, key(kDaGrammarRootHeads)),
                      Item("first_occs", t_self.first_occs_, key(kDaFirstOccs)),
                      Item("last_occs", t_self.last_occs_, key(kDaLastOccs)),
                      NDoc(t_self.n_doc_),
                      Part("get_doc", t_self.backend_),
                      Part("oracle", t_self.oracle_, "oracle_")};
  }

  const SLP* slp_ = nullptr;
  const sdsl::int_vector<>* roots_ = nullptr;
  const RootHeads* root_heads_ = nullptr;
  const typename RootHeads::rank_1_type* rank_ = nullptr;
  const typename RootHeads::select_1_type* select_ = nullptr;
  const OccsBVs* first_occs_ = nullptr;
  const OccsBVs* last_occs_ = nullptr;
  std::size_t n_doc_ = 0;
  TBackend backend_;
  TOracle oracle_;
};

// The DA grammar (dret's, slp/da_grammar.h), its packed rules, where each top-level symbol starts, the
// occurrence bitvectors of every rule; the backend's and the oracle's structures.
template <typename TStorage, uint8_t t_width, typename TBackend, typename TOracle>
void construct(GcdaOccsCore<TStorage, t_width, TBackend, TOracle>& t_core, Config& t_config) {
  using Core = GcdaOccsCore<TStorage, t_width, TBackend, TOracle>;
  using namespace conf;
  EnsureCollection<t_width>(t_config);

  const auto key_slp = KeyName(t_config.keys, kDaGrammarPacked);
  const auto key_heads = KeyName(t_config.keys, kDaGrammarRootHeads);
  const auto key_first = KeyName(t_config.keys, kDaFirstOccs);
  const auto key_last = KeyName(t_config.keys, kDaLastOccs);
  if (!sdsl::cache_file_exists<typename Core::SLP>(key_slp, t_config) ||
      !sdsl::cache_file_exists<typename Core::RootHeads>(key_heads, t_config) ||
      !sdsl::cache_file_exists<typename Core::OccsBVs>(key_first, t_config) ||
      !sdsl::cache_file_exists<typename Core::OccsBVs>(key_last, t_config)) {
    auto event = sdsl::memory_monitor::event(key_first);
    grammar::SLP<> slp;
    std::vector<std::size_t> seq;
    LoadOrBuildDaSlp(t_config, sdsl::cache_file_name<std::vector<int>>(t_config.keys[kDA].get<std::string>(), t_config),
                     slp, seq);

    auto bit_compress = [](sdsl::int_vector<>& v) { sdsl::util::bit_compress(v); };
    typename Core::SLP packed(slp, bit_compress, bit_compress);
    sdsl::store_to_cache(packed, key_slp, t_config, true);

    std::size_t n = 0;
    for (auto v : seq) n += slp.SpanLength(v);
    sdsl::bit_vector heads(n, 0);
    for (std::size_t i = 0, pos = 0; i < seq.size(); pos += slp.SpanLength(seq[i]), ++i) heads[pos] = 1;
    typename Core::RootHeads root_heads(heads);
    sdsl::store_to_cache(root_heads, key_heads, t_config, true);
    typename Core::RootHeads::rank_1_type rank(&root_heads);
    sdsl::store_to_cache(rank, key_heads, t_config, true);
    typename Core::RootHeads::select_1_type select(&root_heads);
    sdsl::store_to_cache(select, key_heads, t_config, true);

    typename Core::OccsBVs first(slp.Variables() + 1), last(slp.Variables() + 1);
    // Every terminal too: a document that occurs only in the top-level sequence is in no rule, and a
    // range of it alone is covered by its terminal, whose bitvectors no rule would have built.
    for (std::size_t i = 0; i <= slp.Sigma(); ++i) BuildFirstAndLastOccs(slp, i, first, last);
    for (auto i = slp.Start(); slp.Sigma() < i; --i) BuildFirstAndLastOccs(slp, i, first, last);
    sdsl::store_to_cache(first, key_first, t_config, true);
    sdsl::store_to_cache(last, key_last, t_config, true);
  }
  EnsureNDoc<t_width>(t_config);
  {
    auto backend = t_core.backend();
    construct(backend, t_config);
  }
  {
    auto oracle = t_core.oracle();
    construct(oracle, t_config);
  }
}

// Every document of the range with its frequency (Weighted).
template <typename TStorage, uint8_t t_width, typename TBackend, typename TOracle>
struct ComputesFrequencies<GcdaOccsCore<TStorage, t_width, TBackend, TOracle>> : std::true_type {};

// A GCDA frequency index: the r-index for the range, then GcdaOccsCore.
template <typename TStorage, typename TAlphabet, typename TCountIdx, typename TBackend, typename TOracle>
using DocFreqIdxGCDAOccs =
    DocListIdx<TStorage, TAlphabet, TCountIdx, GcdaOccsCore<TStorage, TAlphabet::int_width, TBackend, TOracle>>;

}  // namespace dret::freq
