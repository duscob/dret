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

// The length of the int_vector stored in t_file (the BWT: the text length), read from its header.
//
// Not through sdsl::int_vector_buffer, which rewrites the file whenever it is closed: its constructor
// opens the output stream whatever the requested mode, and close() writes the header back. Asking a
// read-only-looking buffer for its size touched the file, harmless in content but it dirtied mtimes,
// which is what "did this stage build anything?" checks read. (sri's constructRIndex did the same
// for every index over an r-index, built or not, until it read the length only when it needed it.)
// Mirrors int_vector_buffer's arithmetic: read_header yields the length in bits and the stored width.
inline std::size_t ReadVectorLength(const std::string& t_file) {
  std::ifstream in(t_file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + t_file + "'");
  uint64_t bits = 0;
  uint8_t width = 0;
  sdsl::int_vector<0>::read_header(bits, width, in);
  if (width == 0) throw std::runtime_error("'" + t_file + "' reports width 0");
  return bits / width;
}

template <uint8_t t_width = 8>
class Factory {
 public:
  struct Index {
    std::shared_ptr<dret::DocListIndex<>> idx;
    std::size_t size = 0;  // serialized
  };

  explicit Factory(dret::Config t_config) : config_{std::move(t_config)} {
    seq_size_ = ReadVectorLength(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_));
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
