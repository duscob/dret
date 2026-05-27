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
// Config.keys. serialize() persists base_ so istream reload can reopen them.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <rlcsa.h>

#include "dret/config.h"
#include "dret/construct_base.h"
#include "dret/index_base.h"
#include "dret/size_report.h"

namespace dret::rmq {

template <typename TStorage = GenericStorage, uint8_t t_width = 8>
class GetDocRLCSA : public IndexBaseWithExternalStorage<TStorage, t_width> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage, t_width>;
  using typename Base::size_type;

  explicit GetDocRLCSA(const TStorage& t_storage) : Base(t_storage) {}
  GetDocRLCSA() = default;

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

  // Serialize the RLCSA sidecar base path so the istream reload path can
  // reopen the sidecar without needing a Config object.
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v,
                      const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    uint64_t len = base_.size();
    out.write(reinterpret_cast<const char*>(&len), sizeof(len));
    out.write(base_.data(), static_cast<std::streamsize>(len));
    return sizeof(len) + len;
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
      auto& in = std::get<std::reference_wrapper<std::istream>>(t_source).get();
      uint64_t len = 0;
      in.read(reinterpret_cast<char*>(&len), sizeof(len));
      base_.resize(len);
      in.read(base_.data(), static_cast<std::streamsize>(len));
    }
    rlcsa_ = std::make_shared<CSA::RLCSA>(base_, /*print=*/false);
    if (!rlcsa_->isOk())
      throw std::runtime_error("GetDocRLCSA: failed to load RLCSA from " + base_);
  }

 private:
  std::shared_ptr<CSA::RLCSA> rlcsa_;
  std::string base_;
};

// Free-function construct(). Builds (or short-circuits when warm) the RLCSA
// sidecar files from the kText cache, remapping dret's \x01 doc-delimiters
// back to the \x00 that RLCSA expects.
template <typename TStorage, uint8_t t_width>
void construct(GetDocRLCSA<TStorage, t_width>&, Config& t_config) {
  using namespace dret::conf;

  if (!cache_file_exists(t_config.keys[kText].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    ConstructText<t_width>(t_config);
  }

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
    throw std::runtime_error("construct(GetDocRLCSA): RLCSA construction failed");
  rlcsa.writeTo(base);
}

}  // namespace dret::rmq
