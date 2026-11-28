//
// RLCSA-backed document-array lookup policy — the paper-faithful baseline for
// PDL (Cobas/Mäkinen/Rossi SPIRE 2020). Uses CSA::RLCSA::locate(range) for
// fast batched SA lookup and getSequenceForPosition for doc-id recovery.
//
// IMPORTANT — PDL ONLY (not RMQ-compatible):
//   RLCSA builds a Generalized-SA where each \x00 sequence separator gets a
//   UNIQUE sort value (rlcsa simpleSuffixSort: `value += zeros; zeros++`),
//   while dret's SA treats all \x01 separators as identical and compares
//   across them. The two SAs cover the same user-suffix SET but in different
//   orders, so dret's SA index k does NOT correspond to RLCSA compact index
//   k-1-D suffix-for-suffix.
//
//   This is fine for PDL (it processes the whole range as a SET, and the SET
//   of doc IDs is preserved between the two orderings), but BREAKS SADA/ILCP
//   cores whose RMQ-driven recursion picks specific positions in dret's SA
//   order. CILCP accidentally works because it doesn't early-stop on marked
//   docs. See docs/pdl_rlcsa_baseline_plan.md for the analysis. The RMQ
//   factory (factories/rmq.h) falls RLCSA through to DA for this reason.
//
// Cache strategy: RLCSA writes its own sidecar files (.rlcsa.array,
// .rlcsa.parameters, .rlcsa.sa_samples) via writeTo(base). These do not
// participate in SDSL typed-cache; the key kRLCSA maps to a file basename via
// Config.keys. serialize() writes the contents of the files the RLCSA reads, so
// a stream holds the whole structure and an istream load needs nothing else.
//

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <rlcsa.h>

#include "dret/archive.h"
#include "dret/collection.h"
#include "dret/backend/backend.h"
#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret::backend {

template <typename TStorage = GenericStorage, uint8_t t_width = 8>
class RLCSA : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit RLCSA(const TStorage& t_storage) : Base(t_storage) {}
  RLCSA() = default;

  // Single-position lookup.
  std::size_t operator()(std::size_t t_i) const {
    const CSA::usint d = rlcsa_->getNumberOfSequences();
    return rlcsa_->getSequenceForPosition(rlcsa_->locate(t_i - 1 - d));
  }

  // Batched range [b, e). Apply the same -(1+D) offset to both bounds.
  template <typename TReport>
  void operator()(std::size_t t_b, std::size_t t_e, TReport& t_report) const {
    if (t_b >= t_e) return;
    const CSA::usint d = rlcsa_->getNumberOfSequences();
    CSA::pair_type range(t_b - 1 - d, t_e - 2 - d);
    CSA::usint* res = rlcsa_->locate(range);
    rlcsa_->getSequenceForPosition(res, CSA::length(range));
    for (CSA::usint i = range.first; i <= range.second; ++i) {
      t_report(static_cast<std::size_t>(res[i - range.first]));
    }
    delete[] res;
  }

  // In-order expansion of [b, e) that stops once f returns false. RMQ x RLCSA
  // is never built (the factory falls back to DA, see above), so this only has
  // to be correct: single lookups, left to right.
  template <typename F>
  void ExpandUntil(std::size_t t_b, std::size_t t_e, F&& f) const {
    for (std::size_t i = t_b; i < t_e; ++i)
      if (!f((*this)(i)))
        return;
  }

  // Suffix-array range of t_pattern in dret's positions, [sp, ep), found by the RLCSA itself: an index
  // over this backend counts with it (as the original PDL does) instead of loading an r-index too.
  // Inverse of the lookups above (dret position i = RLCSA position i - 1 - D); the two suffix arrays
  // order suffixes differently only around the separators, so a pattern's range has the same bounds.
  template <typename TPattern>
  std::pair<std::size_t, std::size_t> Count(const TPattern& t_pattern) const {
    const CSA::pair_type range = rlcsa_->count(std::string(t_pattern.begin(), t_pattern.end()));
    if (CSA::isEmpty(range)) return {0, 0};
    const CSA::usint d = rlcsa_->getNumberOfSequences();
    return {range.first + 1 + d, range.second + 2 + d};
  }

  // The RLCSA reads and writes only files, so serialize() writes the contents of each file it reads,
  // in kFiles order. The serialized size is then the RLCSA's own, and no path enters the stream: it
  // used to hold the files' path, so its bytes depended on the directory the index lived in.
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v,
                      const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written = 0;
    for (const auto& ext : kFiles) written += archive::WriteValue(ReadFile(base_ + ext), out, child, ext);
    sdsl::structure_tree::add_size(child, written);
    return written;
  }

  SizeReport GetSizeReport() const {
    SizeReport r;
    if (rlcsa_)
      append(r, "rlcsa", rlcsa_->reportSize());
    return r;
  }

 protected:
  using typename Base::TSource;

  void loadInner(TSource& t_source, const JSON& t_keys) override {
    if (std::holds_alternative<std::reference_wrapper<Config>>(t_source)) {
      auto& config = std::get<std::reference_wrapper<Config>>(t_source).get();
      base_ = sdsl::cache_file_name(
          t_keys[dret::conf::kRLCSA].template get<std::string>(), config);
    } else {
      // Write the files back to a directory of this index's own, kept while the index lives, so a
      // later serialize() can read them again.
      auto& in = std::get<std::reference_wrapper<std::istream>>(t_source).get();
      files_dir_ = std::make_shared<TempDir>();
      base_ = (files_dir_->path / "rlcsa").string();
      for (const auto& ext : kFiles) {
        std::string bytes;
        archive::ReadValue(bytes, in);
        std::ofstream(base_ + ext, std::ios::binary).write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
      }
    }
    rlcsa_ = std::make_shared<CSA::RLCSA>(base_, /*print=*/false);
    if (!rlcsa_->isOk())
      throw std::runtime_error("RLCSA: failed to load RLCSA from " + base_);
  }

 private:
  // The files CSA::RLCSA(base) reads (rlcsa.cpp), in the order serialize() writes them.
  static constexpr std::array<const char*, 3> kFiles{".rlcsa.array", ".rlcsa.parameters", ".rlcsa.sa_samples"};

  static std::string ReadFile(const std::string& t_path) {
    std::ifstream in(t_path, std::ios::binary);
    if (!in) throw std::runtime_error("RLCSA: cannot read " + t_path);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
  }

  // A directory removed with its files when the last index sharing it goes.
  struct TempDir {
    std::filesystem::path path;
    TempDir() {
      std::string tmpl = (std::filesystem::temp_directory_path() / "dret_rlcsa_XXXXXX").string();
      if (mkdtemp(tmpl.data()) == nullptr) throw std::runtime_error("RLCSA: mkdtemp failed");
      path = tmpl;
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    ~TempDir() {
      std::error_code ec;
      std::filesystem::remove_all(path, ec);
    }
  };

  std::shared_ptr<CSA::RLCSA> rlcsa_;
  std::string base_;
  std::shared_ptr<TempDir> files_dir_;  // set by an istream load: the directory base_ is in
};

// Free-function construct(). Builds (or short-circuits when warm) the RLCSA
// sidecar files from the kText cache, remapping dret's \x01 doc-delimiters
// back to the \x00 that RLCSA expects.
template <typename TStorage, uint8_t t_width>
void construct(RLCSA<TStorage, t_width>&, Config& t_config) {
  using namespace dret::conf;

  EnsureText<t_width>(t_config);

  const auto base = sdsl::cache_file_name(
      t_config.keys[kRLCSA].get<std::string>(), t_config);

  // Short-circuit if sidecar already on disk.
  {
    CSA::RLCSA probe(base, /*print=*/false);
    if (probe.isOk()) return;
  }

  // Load dret's kText cache (int_vector<8>) and remap \x01 → \x00.
  // dret text: doc1\x01doc2\x01...docN\x01\x00  (last byte is the SDSL terminal)
  // RLCSA wants: doc1\x00doc2\x00...docN\x00    (each \x00 marks end-of-sequence)
  // → drop the SDSL terminal and remap the delimiters.
  sdsl::int_vector<8> text;
  sdsl::load_from_cache(text, t_config.keys[kText].get<std::string>(), t_config);

  const std::size_t n = text.size() - 1;  // drop the trailing SDSL \x00 terminal
  std::vector<CSA::uchar> buf(n);
  for (std::size_t i = 0; i < n; ++i) {
    buf[i] = (text[i] == 1) ? CSA::uchar(0) : static_cast<CSA::uchar>(text[i]);
  }
  // buf now ends with \x00 (the remapped last \x01 doc-separator)

  CSA::RLCSA rlcsa(buf.data(), static_cast<CSA::usint>(n), /*block_size=*/32,
                   /*sa_sample_rate=*/128, /*threads=*/1, /*delete_data=*/false);
  if (!rlcsa.isOk())
    throw std::runtime_error("construct(RLCSA): RLCSA construction failed");
  rlcsa.writeTo(base);
}

}  // namespace dret::backend
