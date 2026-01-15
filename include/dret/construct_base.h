//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/17/25.
//

#pragma once

#include <string>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "config.h"

#ifndef REPAIR_EXE
#define REPAIR_EXE nullptr
#endif

namespace dret {

namespace conf {

const std::string KEY_DOC_END = "doc_end";

const std::string KEY_DA = "da";
const std::string KEY_DA_RAW = KEY_DA + "_raw";

}  // namespace conf

//~~~~~~~


template <uint8_t t_width>
void ConstructText(Config& t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructText: width must be `0` for integer alphabet and `8` for byte alphabet");
  sdsl::int_vector<t_width> text;
  auto num_bytes = t_config.data_width / 8;
  load_vector_from_file(text, t_config.data_path, num_bytes);

  if (t_config.data_delim != 0 && std::find(text.begin(), text.end(), 0) != text.end()) {
    throw std::logic_error(std::string("Data file `") + t_config.data_path.string() + "` contains inner zero symbol.");
  }

  if (t_config.data_delim != 1) {
    if (std::find(text.begin(), text.end(), 1) != text.end()) {
      throw std::logic_error(std::string("Data file `") + t_config.data_path.string()
                             + "` contains document delimiter used internally.");
    }

    std::replace(text.begin(), text.end(), t_config.data_delim, static_cast<uint64_t>(1u));
  }

  if (auto last = text.end() - 1; *last == 1) {
    // Data ends with a document delimiter, but the data end symbol (zero) is missing
    text.resize(text.size() + 1);
    text[text.size() - 1] = 0;
  } else {
    if (*last == 0) {
      if (*(last - 1) != 1) {
        // Data ends with the data end symbol (zero), but the previous document delimiter is missing
        text.resize(text.size() + 1);
        text[text.size() - 2] = 1;
        text[text.size() - 1] = 0;
      }
    } else {
      // A document delimiter and the data end symbol (zero) are missing
      text.resize(text.size() + 2);
      text[text.size() - 2] = 1;
      text[text.size() - 1] = 0;
    }
  }

  sdsl::store_to_cache(text, t_config.keys[conf::kText].get<std::string>(), t_config);
}

//~~~~~~~


template <typename II, typename DocBorder, typename DocDelim>
void ConstructDocBorder(II begin, II end, DocBorder& doc_border, const DocDelim& doc_delim, std::size_t size = 0) {
  if (size == 0)
    size = std::distance(begin, end);

  DocBorder tmp_doc_border(size, 0);

  std::size_t i = 0;
  for (auto it = begin; it != end; ++it, ++i) {
    if (*it == doc_delim) {
      tmp_doc_border[i] = 1;
    }
  }

  doc_border.swap(tmp_doc_border);
};

//~~~~~~~


template <uint8_t t_width, typename DocBorder, typename DocDelim>
void ConstructDocBorder(const std::string& data_file, DocBorder& doc_border, const DocDelim& doc_delim) {
  sdsl::int_vector_buffer<t_width> data_buf(data_file, std::ios::in, 1024 * 1024, t_width, true);

  DocBorder tmp_doc_border(data_buf.size(), 0);

  for (std::size_t i = 0; i < data_buf.size(); ++i) {
    if (data_buf[i] == doc_delim) {
      tmp_doc_border[i] = 1;
    }
  }

  doc_border.swap(tmp_doc_border);

  // Previous solution seems faster due to iterator comparison.
  //  ConstructDocBorder(data_buf.begin(), data_buf.end(), doc_border, doc_delim);
}

//~~~~~~~


template <uint8_t t_width = 8, typename BitVector = sdsl::sd_vector<>>
void ConstructDocEnd(Config& t_config, uint8_t kDocDelimiter = 1) {
  static_assert(t_width == 0 or t_width == 8,
                "constructDocEnd: width must be `0` for integer alphabet and `8` for byte alphabet");

  sdsl::bit_vector tmp_doc_endings;
  auto key_text = t_config.keys[conf::kText].get<std::string>();
  ConstructDocBorder<t_width>(sdsl::cache_file_name(key_text, t_config), tmp_doc_endings, kDocDelimiter);

  BitVector doc_endings(tmp_doc_endings);
  auto key_doc_end = t_config.keys[conf::kDocEnds].get<std::string>();
  sdsl::store_to_cache(doc_endings, key_doc_end, t_config);
  sdsl::store_to_cache(doc_endings, key_doc_end, t_config, true);

  typename BitVector::rank_1_type doc_endings_rank(&doc_endings);
  sdsl::store_to_cache(doc_endings_rank, key_doc_end, t_config, true);

  typename BitVector::select_1_type doc_endings_select(&doc_endings);
  sdsl::store_to_cache(doc_endings_select, key_doc_end, t_config, true);
}

//~~~~~~~


template <typename BitVector = sdsl::sd_vector<>>
void ConstructDocArray(sdsl::cache_config& t_config) {
  sdsl::int_vector<> da;

  {
    sdsl::int_vector_buffer<> sa_buf(sdsl::cache_file_name(sdsl::conf::KEY_SA, t_config), std::ios::in, 1024 * 1024);

    BitVector doc_endings;
    load_from_cache(doc_endings, conf::KEY_DOC_END, t_config);
    typename BitVector::rank_1_type doc_endings_rank(&doc_endings);
    auto doc_cnt = doc_endings_rank(doc_endings.size());


    da = sdsl::int_vector<>(sa_buf.size(), 0, sdsl::bits::hi(doc_cnt) + 1);
    for (size_t i = 0; i < sa_buf.size(); ++i) {
      da[i] = doc_endings_rank(sa_buf[i]);
    }

    store_to_cache(da, conf::KEY_DA, t_config);
  }

  {
    std::vector<int> da_raw;
    da_raw.reserve(da.size());
    for (auto&& i : da) {
      da_raw.emplace_back(i);
    }

    auto filepath = cache_file_name(conf::KEY_DA_RAW, t_config);
    sdsl::osfstream out(filepath, std::ios::binary | std::ios::trunc | std::ios::out);
    serialize_vector(da_raw, out);
  }
}

//~~~~~~~


}  // namespace dret
