//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/26.
//
// The grammar of the document array: the one place it is compressed, and its two readings.
//
// EnsureDaGrammarFiles runs the grammar compressor over the DA file (std::vector<int>, built by
// EnsureDA) unless its output is there, and checks the output. The compressor is RePair (irepair,
// an external process; see dret/repair.h for how its binary and memory are chosen), and the output
// is its .R (rules) and .C (top-level sequence) next to the DA file. Every grammar dret builds over
// the DA starts from that output, read two ways:
//   LoadOrBuildDaCnfGrammar  the grammar in CNF (the top-level sequence folded into binary rules),
//                            cached as da_cnf_grammar: the sampled tree, the compact encodings, the
//                            bare grammar of Brute-GCDA and GCDA-differential's tree start from it;
//   LoadOrBuildDaSlp         the rules plus the top-level sequence, cached as da_grammar and
//                            da_grammar_seq: the sampled-cached representation (LightSLP) uses it.
//
// Another compressor plugs in here: a case in EnsureDaGrammarFiles that writes its grammar in the
// same .R/.C form (binary rules; run-length and multi-symbol rules binarized), chosen by a setting in
// dret::Config, plus a term for that setting in the key of every component built from the grammar
// (KeyName / CellKey in slp/slp_components.h), so the two grammars never share a cache file. The
// differential grammar (slp/differential_*.h) is compressed separately, in process, from the
// differences of the DA.
//

#pragma once

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

// For decoding the std::system() wait status below. This pipeline is Linux-only
// (irepair ships as a gcc Makefile), so the POSIX macros are always available.
#include <sys/wait.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/memory_management.hpp>

#include <grammar/re_pair.h>
#include <grammar/slp.h>
#include <grammar/slp_helper.h>

#include "dret/config.h"
#include "dret/repair.h"
#include "dret/slp/slp_components.h"

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

// The compressor's output for the DA file t_datafile, built unless there: RePair's .R and .C next to
// it, registered in the config's file map and checked. The one place the DA's grammar is compressed.
inline void EnsureDaGrammarFiles(Config& t_config, const std::string& t_datafile) {
  if (!std::filesystem::exists(t_datafile + ".R") && repair::kAvailable) {
    const auto filename = std::filesystem::path(t_datafile).filename().string();
    auto event = sdsl::memory_monitor::event("RePair-" + filename);
    RunRePair(t_datafile, t_config.repair);
    t_config.file_map[filename + ".R"] = t_datafile + ".R";
    t_config.file_map[filename + ".C"] = t_datafile + ".C";
  }
  CheckRePairGrammar(t_datafile);
}

//~~~~~~~

// The CNF grammar of the DA, from irepair's output for the DA file (running
// irepair if needed). Cached once per collection as the da_cnf_grammar
// component in grammar::SLP<>'s encoding, the one every builder starts from:
// every GCDA representation, GCDA-differential's sampled tree and GCDA-nolists.
//
// It is read from irepair's output, never re-encoded in process: the in-process
// grammar::RePairBasicEncoder degrades pathologically on some document arrays
// (on concat_1000_001 it ran >45 h without finishing, 80% of samples in
// grammar::searchHash), where irepair took 479 s. See docs/bug_dgcda_repair_hang.md.
inline grammar::SLP<> LoadOrBuildDaCnfGrammar(Config& t_config, const std::string& t_datafile) {
  grammar::SLP<> slp;
  const auto key = KeyName(t_config.keys, conf::kDaCnfGrammar);
  if (sdsl::load_from_cache(slp, key, t_config, true))
    return slp;

  EnsureDaGrammarFiles(t_config, t_datafile);

  auto event = sdsl::memory_monitor::event(key);
  {
    grammar::RePairReader<true> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(slp);
    re_pair_reader.Read(t_datafile, slp_wrapper);
  }
  sdsl::store_to_cache(slp, key, t_config, true);
  return slp;
}

//~~~~~~~


// The SLP over the document array, plus the compact sequence, as parsed from
// irepair's .R/.C. Cached once per collection under bs/sf-independent keys.
//
// irepair itself already ran once per collection -- every call site guards it
// on std::filesystem::exists(t_datafile + ".R"), a path with no (bs,sf) in it.
// What used to repeat was this *parse* of its output: it lived inside the
// per-cell guard in construct() below, so all 20 cells of a grid sweep redid
// it. Measured across the 2026-08 campaign that cost ~71 h, ~8.4% of stage A,
// and it is why revision-mid spent 12.85 h re-deriving one grammar 17 times.
//
// compact_seq round-trips through a bit-compressed int_vector; the values are
// grammar variable ids, so they are non-negative and bounded by the rule count.
inline void LoadOrBuildDaSlp(Config& t_config,
                             const std::string& t_datafile,
                             grammar::SLP<>& t_slp,
                             std::vector<std::size_t>& t_compact_seq) {
  using namespace conf;
  const auto key_slp = KeyName(t_config.keys, kDaGrammar);
  const auto key_seq = KeyName(t_config.keys, kDaGrammarSeq);

  if (sdsl::cache_file_exists<grammar::SLP<>>(key_slp, t_config)
      && sdsl::cache_file_exists<sdsl::int_vector<>>(key_seq, t_config)) {
    sdsl::load_from_cache(t_slp, key_slp, t_config, true);
    sdsl::int_vector<> seq;
    sdsl::load_from_cache(seq, key_seq, t_config, true);
    t_compact_seq.assign(seq.begin(), seq.end());
    return;
  }

  EnsureDaGrammarFiles(t_config, t_datafile);
  auto event = sdsl::memory_monitor::event(key_slp);
  {
    grammar::RePairReader<false> re_pair_reader;
    auto slp_wrapper = grammar::BuildSLPWrapper(t_slp);
    auto report_compact_seq = [&t_compact_seq](const auto& _var) {
      t_compact_seq.emplace_back(_var);
    };
    re_pair_reader.Read(t_datafile, slp_wrapper, report_compact_seq);
  }

  sdsl::store_to_cache(t_slp, key_slp, t_config, true);
  sdsl::int_vector<> seq(t_compact_seq.size());
  std::copy(t_compact_seq.begin(), t_compact_seq.end(), seq.begin());
  sdsl::util::bit_compress(seq);
  sdsl::store_to_cache(seq, key_seq, t_config, true);
}

//~~~~~~~

}  // namespace dret
