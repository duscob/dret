//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// The collection-level structures every index starts from: the text (with dret's document delimiter),
// its suffix array, the document ends, the document array, the number of documents and the LCP array.
// Each Ensure* builds its structure into the cache unless it is there already, after the structures it
// is computed from, so a construct() asks only for what it reads. They are shared by every index over
// the collection and keyed by the common keys alone (no parameter of any index).
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include <sdsl/construct_isa.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/memory_management.hpp>
#include <sdsl/sd_vector.hpp>

#include "dret/config.h"
#include "dret/construct_base.h"

namespace dret {

template <uint8_t t_width>
void EnsureText(Config& t_config) {
  if (!cache_file_exists(t_config.keys[conf::kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<t_width>(t_config);
  }
}

template <uint8_t t_width>
void EnsureSA(Config& t_config) {
  EnsureText<t_width>(t_config);
  if (const auto key = t_config.keys[conf::kSA].get<std::string>(); !cache_file_exists(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    sdsl::construct_sa<t_width>(t_config);
  }
}

// ConstructDocEnd stores the bitvector with and without its type hash; the check is by type, so a
// bitvector of another type is built rather than taken for this one.
template <uint8_t t_width, typename TBvDocEnds = sdsl::sd_vector<>>
void EnsureDocEnds(Config& t_config) {
  EnsureText<t_width>(t_config);
  if (const auto key = t_config.keys[conf::kDocEnds].get<std::string>();
      !sdsl::cache_file_exists<TBvDocEnds>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocEnd<t_width, TBvDocEnds>(t_config);
  }
}

// ConstructDocArray stores the DA only with its type hash, so the check must be by type too, or the DA
// is rebuilt on every call.
template <uint8_t t_width, typename TBvDocEnds = sdsl::sd_vector<>>
void EnsureDA(Config& t_config) {
  EnsureSA<t_width>(t_config);
  EnsureDocEnds<t_width, TBvDocEnds>(t_config);
  if (const auto key = t_config.keys[conf::kDA].get<std::string>();
      !sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    auto event = sdsl::memory_monitor::event(key);
    ConstructDocArray<TBvDocEnds>(t_config);
  }
}

// The number of documents: the document delimiters of the cached text, cached at conf::kNDoc as a
// one-element int_vector<64> so an index can load it without the text or the document ends.
template <uint8_t t_width>
void EnsureNDoc(Config& t_config) {
  EnsureText<t_width>(t_config);
  const auto key = t_config.keys[conf::kNDoc].get<std::string>();
  if (sdsl::cache_file_exists<sdsl::int_vector<64>>(key, t_config)) return;

  sdsl::int_vector_buffer<8> text_buf(sdsl::cache_file_name(t_config.keys[conf::kText].get<std::string>(), t_config));
  std::size_t n_doc = 0;
  for (std::size_t i = 0; i < text_buf.size(); ++i) {
    if (text_buf[i] == 1) ++n_doc;
  }
  sdsl::int_vector<64> n_doc_vec(1, n_doc);
  sdsl::store_to_cache(n_doc_vec, key, t_config, true);
}

// The cached number of documents, built first if missing.
template <uint8_t t_width>
std::size_t ReadNDoc(Config& t_config) {
  EnsureNDoc<t_width>(t_config);
  sdsl::int_vector<64> n_doc_vec;
  sdsl::load_from_cache(n_doc_vec, t_config.keys[conf::kNDoc].get<std::string>(), t_config, true);
  return static_cast<std::size_t>(n_doc_vec[0]);
}

// The LCP array (Kasai, from the SA and its inverse).
template <uint8_t t_width>
void EnsureLCP(Config& t_config) {
  EnsureSA<t_width>(t_config);
  if (!cache_file_exists(sdsl::conf::KEY_LCP, t_config)) {
    auto event = sdsl::memory_monitor::event(std::string(sdsl::conf::KEY_LCP));
    if (!cache_file_exists(sdsl::conf::KEY_ISA, t_config)) {
      sdsl::construct_isa(t_config);
    }
    sdsl::construct_lcp_kasai<t_width>(t_config);
  }
}

// What the indexes over the document array start from: text, SA, document ends, DA and the number of
// documents.
template <uint8_t t_width, typename TBvDocEnds = sdsl::sd_vector<>>
void EnsureCollection(Config& t_config) {
  EnsureDA<t_width, TBvDocEnds>(t_config);
  EnsureNDoc<t_width>(t_config);
}

}  // namespace dret
