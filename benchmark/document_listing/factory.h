//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/3/21.
//
// The query-mode factory: the storage the indexes of a query run share, the text length, and each
// cell's index (factories/cells.h), built if its cache files are missing and loaded on first use.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include "sr-index/index_base.h"

#include "dret/config.h"
#include "dret/doc_list/doc_list_base.h"

#include "factories/cells.h"

template <uint8_t t_width = 8>
class Factory {
 public:
  struct Index {
    std::shared_ptr<dret::DocListIndex<>> idx;
    std::size_t size = 0;  // serialized
  };

  explicit Factory(dret::Config t_config) : config_{std::move(t_config)} {
    // Sequence length, read straight from the BWT's int_vector header.
    //
    // This used to be `int_vector_buffer<t_width> buf(...); seq_size_ = buf.size();`
    // which REWROTE the file on every query run. sdsl::int_vector_buffer's
    // constructor forces the output stream open regardless of the requested
    // mode --- `m_ofile.open(m_filename, mode | std::ios::out | std::ios::binary)`
    // --- and its destructor calls close(), which unconditionally does
    // write_block() and rewrites the header. So merely asking a read-only-looking
    // buffer for its size touched bwt_data.sdsl every time a Factory was built.
    //
    // Harmless in content (the BWT is deterministic, and the 8 collections whose
    // copies happened to be read-only completed identically with the write
    // failing), but it dirtied mtimes, which is what makes "did the query stage
    // build anything?" checks unreliable --- they compare filenames, not mtimes.
    //
    // Mirrors int_vector_buffer's own arithmetic: read_header yields the length
    // in BITS and the stored width, and m_size = bits / width.
    {
      const auto bwt_file = sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_);
      std::ifstream bwt_in(bwt_file, std::ios::binary);
      if (!bwt_in) {
        throw std::runtime_error("Factory: cannot open BWT file '" + bwt_file + "'");
      }
      uint64_t bits = 0;
      uint8_t stored_width = 0;
      sdsl::int_vector<0>::read_header(bits, stored_width, bwt_in);
      if (stored_width == 0) {
        throw std::runtime_error("Factory: BWT '" + bwt_file + "' reports width 0");
      }
      seq_size_ = bits / stored_width;
    }
  }

  // The index of t_cell, made once per run (a benchmark runs once per repetition).
  Index Get(const bench::factories::Cell& t_cell) {
    if (auto it = indexes_.find(t_cell.name); it != indexes_.end()) return it->second;
    auto cell_index = t_cell.make(std::ref(storage_));
    if (t_cell.query_constructs) cell_index.construct(config_);
    cell_index.load(config_);
    return indexes_[t_cell.name] = {cell_index.index, cell_index.serialized_size()};
  }

  [[nodiscard]] std::size_t SequenceSize() const { return seq_size_; }

  // Occurrences of t_pattern in the text: the size of its suffix-array range in the r-index. The
  // r-index is loaded on the first call, not here: it needs the r-index sample / mark cache files,
  // which a run over the other families need not have.
  std::size_t CountOccurrences(const std::string& t_pattern) {
    if (!count_occurrences_) count_occurrences_ = bench::factories::OccurrenceCounter(std::ref(storage_), config_);
    return count_occurrences_(t_pattern);
  }

 private:
  dret::Config config_;
  std::size_t seq_size_ = 0;
  sri::GenericStorage storage_;
  std::map<std::string, Index> indexes_;
  std::function<std::size_t(const std::string&)> count_occurrences_;
};
