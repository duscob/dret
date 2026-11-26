//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Frequency oracles: the frequency of a document in a suffix-array range from the leftmost and the
// rightmost position the document has in it. A pattern's occurrences in a document are contiguous in
// the document's own suffix array, and the suffix array of the collection orders a document's suffixes
// as the document's suffix array does (the delimiter sorts below every symbol), so the two extremes
// bound them. The 2020 frequency code had both:
//
//   WtRankOracle   Rank(r + 1, d) - Rank(l, d), from a backend that counts (backend::WaveletDA)
//   DocIsaOracle   |ISA_d(SA[l]) - ISA_d(SA[r])| + 1, with the per-document inverse suffix arrays in a
//                  differential grammar, by text position, and SA[] from a backend that locates
//                  (backend::DiffSA)
//
// An oracle is a component (it loads, serializes and sizes what it stores) with
//
//   Freq(backend, d, l, r)   the frequency of d, whose leftmost and rightmost positions are l <= r
//
// and a construct(oracle, config).
//

#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include <sdsl/construct_sa.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>

#include "dret/collection.h"
#include "dret/config.h"
#include "dret/fields.h"
#include "dret/index_base.h"
#include "dret/size_report.h"
#include "dret/slp/differential_slp.h"
#include "dret/slp/slp_components.h"

namespace dret::freq {

// From the backend's rank: it stores nothing.
template <typename TStorage = GenericStorage, uint8_t t_width = 8>
class WtRankOracle : public FieldsComponent<WtRankOracle<TStorage, t_width>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<WtRankOracle, TStorage, t_width>;

  explicit WtRankOracle(const TStorage& t_storage) : Base(t_storage) {}

  WtRankOracle() = default;

  template <typename TBackend>
  std::size_t Freq(const TBackend& t_backend, std::size_t t_doc, std::size_t t_l, std::size_t t_r) const {
    return t_backend.Rank(t_r + 1, t_doc) - t_backend.Rank(t_l, t_doc);
  }

 protected:
  friend Base;

  template <typename TSelf>
  static auto Fields(TSelf&) {
    return std::tuple{};
  }
};

template <typename TStorage, uint8_t t_width>
void construct(WtRankOracle<TStorage, t_width>&, Config&) {}

// The per-document inverse suffix arrays, concatenated by text position: at the position of the k-th
// symbol of a document, the rank of that document's k-th suffix among its suffixes (0 at the
// delimiters and the terminator). Stored as a differential grammar with samples every block.
template <typename TStorage = GenericStorage, uint8_t t_width = 8, typename TDSLP = DifferentialSLP<>>
class DocIsaOracle : public FieldsComponent<DocIsaOracle<TStorage, t_width, TDSLP>, TStorage, t_width> {
 public:
  using Base = FieldsComponent<DocIsaOracle, TStorage, t_width>;
  using typename Base::size_type;

  explicit DocIsaOracle(const TStorage& t_storage, uint32_t t_block_size = 256)
      : Base(t_storage), block_size_(t_block_size) {}

  DocIsaOracle() = default;

  template <typename TBackend>
  std::size_t Freq(const TBackend& t_backend, std::size_t /*t_doc*/, std::size_t t_l, std::size_t t_r) const {
    if (t_l == t_r) return 1;
    const auto a = Isa(t_backend.Locate(t_l));
    const auto b = Isa(t_backend.Locate(t_r));
    return (a < b ? b - a : a - b) + 1;
  }

  uint32_t block_size() const { return block_size_; }

  static std::string Key(const JSON& t_keys, uint32_t t_block_size) {
    return PrefixedKey(t_keys, conf::kBlk, KeyName(t_keys, conf::kDocIsa), t_block_size);
  }

 protected:
  friend Base;

  std::size_t Isa(std::size_t t_text_pos) const {
    std::size_t value = 0;
    auto report = [&value](auto v) { value = static_cast<std::size_t>(v); };
    ExpandSLP(*isa_, t_text_pos, t_text_pos + 1, report);
    return value;
  }

  template <typename TSelf>
  static auto Fields(TSelf& t_self) {
    return std::tuple{Item("doc_isa", t_self.isa_, [&t_self](const JSON& k) { return Key(k, t_self.block_size_); })};
  }

  const TDSLP* isa_ = nullptr;
  uint32_t block_size_ = 256;
};

// The text, a document at a time: the suffix array of its symbols (the delimiter is the end, which
// sorts first), inverted into the document's stretch of the array.
template <typename TStorage, uint8_t t_width, typename TDSLP>
void construct(DocIsaOracle<TStorage, t_width, TDSLP>& t_oracle, Config& t_config) {
  EnsureText<t_width>(t_config);
  const auto key = DocIsaOracle<TStorage, t_width, TDSLP>::Key(t_config.keys, t_oracle.block_size());
  if (sdsl::cache_file_exists<TDSLP>(key, t_config)) return;
  auto event = sdsl::memory_monitor::event(key);

  sdsl::int_vector<8> text;
  sdsl::load_from_cache(text, t_config.keys[conf::kText].get<std::string>(), t_config);
  sdsl::int_vector<> isa(text.size(), 0, sdsl::bits::hi(text.size()) + 1);
  std::vector<unsigned char> doc;
  for (std::size_t begin = 0; begin < text.size();) {
    std::size_t end = begin;
    while (end < text.size() && text[end] > 1) ++end;
    if (end == text.size()) break;  // the terminator
    doc.assign(text.begin() + begin, text.begin() + end);
    doc.push_back(0);
    sdsl::int_vector<> sa(end - begin, 0, 64);
    sdsl::algorithm::calculate_sa(doc.data(), end - begin, sa);
    for (std::size_t r = 0; r < sa.size(); ++r) isa[begin + sa[r]] = r;
    begin = end + 1;
  }

  TDSLP dslp;
  dslp.Compute(isa, t_oracle.block_size());
  sdsl::store_to_cache(dslp, key, t_config, true);
}

}  // namespace dret::freq
