//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/17/25.
//

#pragma once

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>

// For decoding the std::system() wait status below. This pipeline is Linux-only
// (irepair ships as a gcc Makefile), so the POSIX macros are always available.
#include <sys/wait.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "dret/config.h"
#include "dret/repair.h"

// The two balanced irepair builds, named for the width of a sequence position:
// REPAIR_EXE_BAL32 is upstream's `bal/irepair`, REPAIR_EXE_BAL64 its
// `large/bal/irepair`. RunRePair chooses between them per document array.
#ifndef REPAIR_EXE_BAL32
#define REPAIR_EXE_BAL32 nullptr
#endif

#ifndef REPAIR_EXE_BAL64
#define REPAIR_EXE_BAL64 nullptr
#endif

namespace dret {

namespace conf {

const std::string KEY_DOC_END = "doc_end";

const std::string KEY_DA = "da";
const std::string KEY_DA_RAW = KEY_DA + "_raw";

}  // namespace conf

//~~~~~~~


// Cache-integrity guards. See docs/bug_empty_da_page_big.md.
//
// A construct that dies while an sdsl::int_vector_buffer is open leaves a file
// whose header still reads zero -- close() is what back-patches the real header
// -- but which otherwise holds gigabytes of entirely plausible data. Nothing in
// the cache layer validates that on load, so a later run reads size 0, builds an
// empty document array, hands it to RePair, and segfaults in the grammar reader
// a long way from the actual fault. These guards make it fail where it breaks.

inline void CheckNonEmptyCacheFile(std::size_t t_size, const std::string& t_path, const std::string& t_what) {
  if (t_size == 0) {
    throw std::runtime_error(
        t_what + " is empty: '" + t_path +
        "'. A cached artifact of size zero usually means an earlier construct was "
        "killed while writing it, leaving real data behind an unpatched header. "
        "Delete it and rebuild -- do not trust the rest of this cache.");
  }
}

// RePair runs as an external process, and its exit status used to be discarded
// at every call site. Check the input, the status, and the outputs here instead.
// Precedent for the rc check: pdl/set_codecs.h, for vnmextract.
namespace repair {

// Whether this build can run RePair at all. The call sites skip grammar
// construction when it cannot, so they ask this rather than testing one
// particular binary -- which was the old `&& REPAIR_EXE` guard, and would now
// silently mean "is the 32-bit one configured".
inline constexpr bool kAvailable = (REPAIR_EXE_BAL32 != nullptr) && (REPAIR_EXE_BAL64 != nullptr);

}  // namespace repair

inline void RunRePair(const std::string& t_datafile, const repair::Options& t_options = {}) {
  std::error_code ec;
  const auto in_size = std::filesystem::file_size(t_datafile, ec);
  if (ec || in_size == 0) {
    throw std::runtime_error(
        "RePair: refusing to compress '" + t_datafile +
        "': " + (ec ? ec.message() : std::string("the file is empty")) +
        ". An empty document array is the signature of a poisoned suffix-array cache.");
  }

  // Which binary, and with what <MB>, is decided from the size of *this* array.
  // The policy itself lives in PlanRePair so it can be tested directly.
  const auto elems = static_cast<std::uintmax_t>(in_size) / sizeof(int);
  const auto plan = repair::PlanRePair(elems, t_options);

  const char* exe = plan.use_64bit ? REPAIR_EXE_BAL64 : REPAIR_EXE_BAL32;
  if (exe == nullptr) {
    throw std::runtime_error(std::string("RePair: the ") + (plan.use_64bit ? "bal64" : "bal32") +
                             " binary is not configured in this build (REPAIR_EXE_" +
                             (plan.use_64bit ? "BAL64" : "BAL32") + " is unset).");
  }

  if (plan.lean_path) {
    // Below the fast-path threshold irepair emits a different grammar. It stays
    // valid and is reproducible, but it will not match the rest of a cache built
    // at the derived value, so this must never pass unremarked.
    std::cerr << "RePair: running with " << *plan.mb << " MB, below the "
              << repair::FastPathMB(elems) << " MB this array needs for the fast path (peak RSS "
              << "would be about " << repair::FastPathPeakMB(elems) << " MB). The memory-lean "
              << "pass produces a different grammar from one built at the derived value: "
              << t_datafile << std::endl;
  }

  std::string cmd = exe + (" " + t_datafile);
  if (plan.mb.has_value()) {
    cmd += " " + std::to_string(*plan.mb);
  }

  const int rc = std::system(cmd.c_str());
  if (rc != 0) {
    // std::system yields a wait status, not an exit code -- report the decoded
    // form, because "killed by signal 9" (the OOM killer, the likely fate of
    // irepair on the largest document arrays) and "exited 1" need different fixes.
    std::string how;
    if (rc == -1) {
      how = "could not be started";
    } else if (WIFSIGNALED(rc)) {
      how = "was killed by signal " + std::to_string(WTERMSIG(rc));
    } else if (WIFEXITED(rc)) {
      how = "exited with code " + std::to_string(WEXITSTATUS(rc));
    } else {
      how = "returned wait status " + std::to_string(rc);
    }
    throw std::runtime_error("RePair: '" + cmd + "' " + how + ".");
  }
}

// Validate the .R/.C pair whether it was just built or picked up from cache --
// a poisoned grammar is otherwise reused forever, since the call sites skip
// RePair whenever .R merely exists.
//
// Size is only decisive for .C: irepair emits a 4-byte .R even for empty input,
// so an empty grammar is recognised by its 0-byte sequence file, not its rules.
inline void CheckRePairGrammar(const std::string& t_datafile) {
  std::error_code ec;
  for (const auto* ext : {".R", ".C"}) {
    const auto path = t_datafile + ext;
    const auto size = std::filesystem::file_size(path, ec);
    if (ec) {
      throw std::runtime_error("RePair: cannot stat '" + path + "': " + ec.message() +
                               ". The grammar for this document array is missing or unreadable.");
    }
    if (size == 0) {
      throw std::runtime_error("RePair: '" + path +
                               "' is empty, so the grammar encodes nothing. Delete the .R/.C pair "
                               "for this document array and rebuild it.");
    }
  }
}

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
void ConstructDocBorder(const std::string& data_file, DocBorder& doc_border, const DocDelim& doc_delim, bool is_plain) {
  sdsl::int_vector_buffer<t_width> data_buf(data_file, std::ios::in, 1024 * 1024, t_width, is_plain);

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
  ConstructDocBorder<t_width>(sdsl::cache_file_name(key_text, t_config), tmp_doc_endings, kDocDelimiter, false);

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
void ConstructDocArray(Config& t_config) {
  sdsl::int_vector<> da;
  auto key_da = t_config.keys[conf::kDA].get<std::string>();

  {
    auto key_sa = t_config.keys[conf::kSA].get<std::string>();
    auto sa_path = sdsl::cache_file_name(key_sa, t_config);
    sdsl::int_vector_buffer<> sa_buf(sa_path, std::ios::in);
    CheckNonEmptyCacheFile(sa_buf.size(), sa_path, "suffix array");

    BitVector doc_endings;
    auto key_doc_end = t_config.keys[conf::kDocEnds].get<std::string>();
    load_from_cache(doc_endings, key_doc_end, t_config, true);
    typename BitVector::rank_1_type doc_endings_rank(&doc_endings);
    auto doc_cnt = doc_endings_rank(doc_endings.size());


    da = sdsl::int_vector<>(sa_buf.size(), 0, sdsl::bits::hi(doc_cnt) + 1);
    for (size_t i = 0; i < sa_buf.size(); ++i) {
      da[i] = doc_endings_rank(sa_buf[i]);
    }

    sdsl::store_to_cache(da, key_da, t_config, true);
  }

  {
    std::vector<int> da_raw;
    da_raw.reserve(da.size());
    for (auto&& i : da) {
      da_raw.emplace_back(i);
    }

    auto filepath = sdsl::cache_file_name<decltype(da_raw)>(key_da, t_config);
    sdsl::osfstream out(filepath, std::ios::binary | std::ios::trunc | std::ios::out);
    serialize_vector(da_raw, out);
    sdsl::register_cache_file<decltype(da_raw)>(key_da, t_config);
  }
}

//~~~~~~~


}  // namespace dret
