//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/17/25.
//

#pragma once

#include <string>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector_buffer.hpp>

namespace dret {

namespace conf {

const std::string KEY_DOC_END = "doc_end";

}  // namespace conf

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


template <uint8_t t_width = 8, typename DocBorder, typename DocDelim>
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
void ConstructDocEnd(sdsl::cache_config& t_config, uint8_t kDocDelimiter = 2) {
  static_assert(t_width == 0 or t_width == 8,
                "constructDocEnd: width must be `0` for integer alphabet and `8` for byte alphabet");

  sdsl::bit_vector tmp_doc_endings;
  ConstructDocBorder<t_width>(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, t_config), tmp_doc_endings, kDocDelimiter);
  BitVector doc_endings(tmp_doc_endings);

  sdsl::store_to_cache(doc_endings, conf::KEY_DOC_END, t_config);
}


}  // namespace dret
