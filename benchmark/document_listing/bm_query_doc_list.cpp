//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/3/21.
//

#include <bitset>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>

#include "../bm_query.h"
#include "../bm_query_set.h"
#include "../bm_size_counters.h"

#include "dret/size_report.h"

#include "factory.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_string(data_name, "data", "Data file basename.");
DEFINE_int32(data_width, 8, "Data width in bits: 8, 16, 32 or 64.");
DEFINE_int32(pattern_delim, '\n', "Pattern delimiter.");

DEFINE_int32(min_s, 4, "Minimum sampling parameter s.");
DEFINE_int32(max_s, 128, "Maximum sampling parameter s.");

DEFINE_int32(min_block_size, 512, "Minimum block size for DocListGCDA (power of 2).");
DEFINE_int32(max_block_size, 512, "Maximum block size for DocListGCDA (power of 2).");
DEFINE_int32(min_storing_factor, 4, "Minimum storing factor for DocListGCDA (power of 2).");
DEFINE_int32(max_storing_factor, 4, "Maximum storing factor for DocListGCDA (power of 2).");
DEFINE_string(rmq_get_doc_variants, "da,slp", "RMQ GetDoc variants to run: comma-separated da,slp,dslp.");

DEFINE_string(gcda_slp_variants,
              "default,compact_bp,compact_louds",
              "GCDA TSLP variants: comma-separated default,compact_bp,compact_louds.");

DEFINE_bool(report_stats, false, "Report statistics for benchmark (mean, median, ...).");
DEFINE_int32(reps, 10, "Repetitions for the locate query benchmark.");
DEFINE_double(min_time, 0, "Minimum time (seconds) for the locate query micro benchmark.");

DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

//~~~~~~~

std::vector<Factory<>::GetDocEnum> ParseGetDocVariants(const std::string& value) {
  std::vector<Factory<>::GetDocEnum> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "da") {
      variants.push_back(Factory<>::GetDocEnum::DA);
    } else if (item == "slp") {
      variants.push_back(Factory<>::GetDocEnum::SLP);
    } else if (item == "dslp") {
      variants.push_back(Factory<>::GetDocEnum::DSLP);
    } else if (!item.empty()) {
      throw std::invalid_argument("Unknown --rmq_get_doc_variants item: " + item);
    }
  }
  if (variants.empty())
    variants.push_back(Factory<>::GetDocEnum::DA);
  return variants;
}

const char* GetDocName(Factory<>::GetDocEnum variant) {
  switch (variant) {
    case Factory<>::GetDocEnum::DA:
      return "DA";
    case Factory<>::GetDocEnum::SLP:
      return "SLP";
    case Factory<>::GetDocEnum::DSLP:
      return "DSLP";
  }
  return "UNKNOWN";
}

std::vector<Factory<>::GCDASLPVariant> ParseGCDASLPVariants(const std::string& value) {
  std::vector<Factory<>::GCDASLPVariant> variants;
  std::stringstream ss(value);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item == "default") {
      variants.push_back(Factory<>::GCDASLPVariant::Default);
    } else if (item == "compact_bp") {
      variants.push_back(Factory<>::GCDASLPVariant::CompactBP);
    } else if (item == "compact_louds") {
      variants.push_back(Factory<>::GCDASLPVariant::CompactLOUDS);
    } else if (!item.empty()) {
      throw std::invalid_argument("Unknown --gcda_slp_variants item: " + item);
    }
  }
  if (variants.empty())
    variants.push_back(Factory<>::GCDASLPVariant::Default);
  return variants;
}

const char* GCDASLPVariantName(Factory<>::GCDASLPVariant variant) {
  switch (variant) {
    case Factory<>::GCDASLPVariant::Default:
      return "Default";
    case Factory<>::GCDASLPVariant::CompactBP:
      return "CompactBP";
    case Factory<>::GCDASLPVariant::CompactLOUDS:
      return "CompactLOUDS";
  }
  return "UNKNOWN";
}

//~~~~~~~


class DocListResult {
 public:
  virtual ~DocListResult() = default;

  virtual void operator()(std::size_t t_doc) = 0;

  virtual void operator()() {}

  virtual void Print(std::ostream& t_os) const = 0;
};

class DocListResultVector : public DocListResult {
 public:
  void operator()(std::size_t t_doc) override {
    result_.emplace_back(t_doc);
  }

  void operator()() override {
    sort(result_.begin(), result_.end());
    result_.erase(unique(result_.begin(), result_.end()), result_.end());
  }

  auto size() const {
    return result_.size();
  }

  void Print(std::ostream& t_os) const override {
    for (const auto& item : result_) {
      t_os << item << '\n';
    }
  }

 private:
  std::vector<std::size_t> result_;
};

void PrintResults(std::ostream& t_os, const DocListResultVector& t_result) {
  t_result.Print(t_os);
}

class DocListResultBitvector : public DocListResult {
 public:
  explicit DocListResultBitvector(std::size_t t_size) : result_(t_size) {}

  void operator()(std::size_t t_doc) override {
    result_[t_doc] = true;
  }

  void Print(std::ostream& t_os) const override {
    for (int i = 0; i < result_.size(); ++i) {
      if (result_[i])
        t_os << i << '\n';
    }
  }

 private:
  sdsl::bit_vector result_;
};

class DocListResultVectorBool : public DocListResult {
 public:
  explicit DocListResultVectorBool(std::size_t t_size) : result_(t_size) {}

  void operator()(std::size_t t_doc) override {
    result_[t_doc] = true;
  }

  void Print(std::ostream& t_os) const override {
    for (int i = 0; i < result_.size(); ++i) {
      if (result_[i])
        t_os << i << '\n';
    }
  }

 private:
  std::vector<bool> result_;
};

// Benchmark Queries on Document Listing Index
auto BM_QueryDocList = [](benchmark::State& t_state,  //
                          auto* t_factory,
                          const auto& t_config,
                          const auto& t_patterns,
                          auto t_seq_size) {
  auto [idx, idx_size] = t_factory->Make(t_config.config);

  for (auto _ : t_state) {
    for (const auto& pattern : t_patterns) {
      auto result = t_config.create_result();
      idx->Search(pattern, std::ref(*result));

      (*result)();

      delete result;
    }
  }

  SetupDefaultCounters(t_state);
  t_state.counters["Size(bytes)"] = idx_size;
  t_state.counters["Bits_x_Symbol"] = idx_size * 8.0 / t_seq_size;
  t_state.counters["Patterns"] = t_patterns.size();
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);

  if (idx) {
    const auto sizes = idx->GetSizeReport();
    appendCounters(t_state, sizes);
    t_state.counters["total_index_bytes"] = static_cast<double>(dret::totalBytes(sizes));
  }
};

auto BM_PrintQueryDocList = [](benchmark::State& t_state,
                               const auto& t_idx_name,
                               auto* t_factory,
                               const auto& t_config,
                               const auto& t_patterns) {
  auto [idx, idx_size] = t_factory->Make(t_config.config);

  for (auto _ : t_state) {
    std::ofstream out(std::string("result-") + t_idx_name + ".txt");
    for (const auto& pattern : t_patterns) {
      out << pattern << std::endl;

      auto result = t_config.create_result();
      idx->Search(pattern, std::ref(*result));

      (*result)();

      result->Print(out);

      delete result;
    }
  }

  SetupDefaultCounters(t_state);
};

int main(int argc, char* argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_patterns.empty() || FLAGS_data_name.empty() || FLAGS_data_dir.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  // Query patterns
  std::vector<std::string> patterns;
  {
    std::ifstream pattern_file(FLAGS_patterns.c_str(), std::ios_base::binary);
    if (!pattern_file) {
      std::cerr << "ERROR: Failed to open patterns file!" << std::endl;
      return 3;
    }

    std::string buf;
    while (std::getline(pattern_file, buf)) {
      if (buf.empty())
        continue;

      patterns.emplace_back(buf);
    }
    pattern_file.close();
  }
  auto patterns_encoded =
      ReadPatterns<dret::Alphabet<>::string_type>(FLAGS_patterns, FLAGS_data_width / 8, FLAGS_pattern_delim);

  // Benchmarks configs
  dret::Config config(FLAGS_data_name, FLAGS_data_dir, sri::SDSL_LIBDIVSUFSORT, true);
  std::vector<Factory<>::GetDocEnum> rmq_get_doc_variants;
  std::vector<Factory<>::GCDASLPVariant> gcda_slp_variants;
  try {
    rmq_get_doc_variants = ParseGetDocVariants(FLAGS_rmq_get_doc_variants);
    gcda_slp_variants = ParseGCDASLPVariants(FLAGS_gcda_slp_variants);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  Factory<> factory(config);

  struct Config {
    const char* name;
    Factory<>::Config config;
    std::function<DocListResult*()> create_result;
  };

  auto create_result_vector = []() {
    return new DocListResultVector();
  };
  auto create_result_bitvector = [n_docs = factory.NDocs()]() {
    return new DocListResultBitvector(n_docs);
  };
  auto create_result_vector_bool = [n_docs = factory.NDocs()]() {
    return new DocListResultVectorBool(n_docs);
  };

  std::vector<Config> index_configs = {
      {"Brute-R-Index-V", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, create_result_vector},
      {"Brute-R-Index-Bv", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, create_result_bitvector},
      {"Brute-R-Index-VB", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, create_result_vector_bool},
      {"Brute-SR-Index-V", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, create_result_vector},
      {"Brute-SR-Index-Bv", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, create_result_bitvector},
      {"Brute-SR-Index-VB", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, create_result_vector_bool},
  };

  auto seq_size = factory.SequenceSize();
  std::string print_bm_prefix = "Print-";
  for (const auto& idx_config : index_configs) {
    benchmark::RegisterBenchmark(idx_config.name, BM_QueryDocList, &factory, idx_config, patterns, seq_size);

    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.name;
      benchmark::RegisterBenchmark(
          print_bm_name.c_str(), BM_PrintQueryDocList, idx_config.name, &factory, idx_config, patterns);
    }
  }

  std::vector<IndexConfig<Factory<>::Config>> idx_configs = {
      {"Brute-RIndex", Factory<>::Config{Factory<>::IndexEnum::BRUTE_R_INDEX}, false},
      {"Brute-SRIndex", Factory<>::Config{Factory<>::IndexEnum::BRUTE_SR_INDEX}, true},
  };

  for (const auto get_doc : rmq_get_doc_variants) {
    if (get_doc == Factory<>::GetDocEnum::DA) {
      idx_configs.push_back({"SADA-DA", Factory<>::Config{Factory<>::IndexEnum::SADA}, false});
      idx_configs.push_back({"ILCP-DA", Factory<>::Config{Factory<>::IndexEnum::ILCP}, false});
      idx_configs.push_back({"CILCP-DA", Factory<>::Config{Factory<>::IndexEnum::CILCP}, false});
    }
  }

  for (int64_t bs = FLAGS_min_block_size; bs <= FLAGS_max_block_size; bs *= 2) {
    for (int64_t sf = FLAGS_min_storing_factor; sf <= FLAGS_max_storing_factor; sf *= 2) {
      // GCDA — one entry per requested TSLP variant (Default, CompactBP, CompactLOUDS).
      for (const auto gcda_slp : gcda_slp_variants) {
        auto suffix = (gcda_slp == Factory<>::GCDASLPVariant::Default)
                          ? std::string{}
                          : std::string("-") + GCDASLPVariantName(gcda_slp);
        auto name = "DocListGCDA" + suffix + "-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
        Factory<>::Config cfg{Factory<>::IndexEnum::GCDA,
                              0,
                              static_cast<uint32_t>(bs),
                              static_cast<float>(sf),
                              Factory<>::GetDocEnum::DA,
                              gcda_slp};
        idx_configs.push_back({name, cfg, false});
      }

      auto dgcda_name = "DocListDGCDA-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
      Factory<>::Config dgcda_cfg{Factory<>::IndexEnum::DGCDA, 0, static_cast<uint32_t>(bs), static_cast<float>(sf)};
      idx_configs.push_back({dgcda_name, dgcda_cfg, false});

      auto dgcda_otf_name = "DocListDGCDA-OTF-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
      Factory<>::Config dgcda_otf_cfg{
          Factory<>::IndexEnum::DGCDA_OTF, 0, static_cast<uint32_t>(bs), static_cast<float>(sf)};
      idx_configs.push_back({dgcda_otf_name, dgcda_otf_cfg, false});

      auto dgcda_crl_name = "DocListDGCDA-CRL-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
      Factory<>::Config dgcda_crl_cfg{
          Factory<>::IndexEnum::DGCDA_CRL, 0, static_cast<uint32_t>(bs), static_cast<float>(sf)};
      idx_configs.push_back({dgcda_crl_name, dgcda_crl_cfg, false});

      auto dgcda_ev_name = "DocListDGCDA-EV-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
      Factory<>::Config dgcda_ev_cfg{
          Factory<>::IndexEnum::DGCDA_EV, 0, static_cast<uint32_t>(bs), static_cast<float>(sf)};
      idx_configs.push_back({dgcda_ev_name, dgcda_ev_cfg, false});

      auto dgcda_dv_name = "DocListDGCDA-DV-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
      Factory<>::Config dgcda_dv_cfg{
          Factory<>::IndexEnum::DGCDA_DV, 0, static_cast<uint32_t>(bs), static_cast<float>(sf)};
      idx_configs.push_back({dgcda_dv_name, dgcda_dv_cfg, false});

      auto dgcda_vv_name = "DocListDGCDA-VV-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
      Factory<>::Config dgcda_vv_cfg{
          Factory<>::IndexEnum::DGCDA_VV, 0, static_cast<uint32_t>(bs), static_cast<float>(sf)};
      idx_configs.push_back({dgcda_vv_name, dgcda_vv_cfg, false});

      for (const auto get_doc : rmq_get_doc_variants) {
        if (get_doc != Factory<>::GetDocEnum::SLP && get_doc != Factory<>::GetDocEnum::DSLP)
          continue;

        // For RMQ-SLP, fan out across the requested gcda_slp variants so the
        // SLP-backed RMQ index reuses the matching GCDA build's SLP cache.
        // For RMQ-DSLP, the GCDA TSLP axis is irrelevant — only register one
        // entry (Default) per (bs, sf, dslp).
        const auto inner_variants =
            (get_doc == Factory<>::GetDocEnum::SLP)
                ? gcda_slp_variants
                : std::vector<Factory<>::GCDASLPVariant>{Factory<>::GCDASLPVariant::Default};

        for (const auto gcda_slp : inner_variants) {
          auto variant_suffix = (get_doc == Factory<>::GetDocEnum::SLP &&
                                 gcda_slp != Factory<>::GCDASLPVariant::Default)
                                    ? std::string("-") + GCDASLPVariantName(gcda_slp)
                                    : std::string{};
          const auto suffix = std::string("-") + GetDocName(get_doc) + variant_suffix +
                              "-bs" + std::to_string(bs) + "-sf" + std::to_string(sf);
          Factory<>::Config sada_cfg{Factory<>::IndexEnum::SADA, 0,
                                     static_cast<uint32_t>(bs), static_cast<float>(sf), get_doc, gcda_slp};
          Factory<>::Config ilcp_cfg{Factory<>::IndexEnum::ILCP, 0,
                                     static_cast<uint32_t>(bs), static_cast<float>(sf), get_doc, gcda_slp};
          Factory<>::Config cilcp_cfg{Factory<>::IndexEnum::CILCP, 0,
                                      static_cast<uint32_t>(bs), static_cast<float>(sf), get_doc, gcda_slp};
          idx_configs.push_back({"SADA" + suffix, sada_cfg, false});
          idx_configs.push_back({"ILCP" + suffix, ilcp_cfg, false});
          idx_configs.push_back({"CILCP" + suffix, cilcp_cfg, false});
        }
      }
    }
  }

  QueryBenchmarkConfig query_bm_config{FLAGS_report_stats, FLAGS_reps, FLAGS_min_time, FLAGS_print_result};

  auto factory_result_vector_int = [&factory](const auto& tt_config) {
    auto idx = factory.MakeIndex(tt_config);

    auto query = [idx = idx.idx](const auto& tt_pattern) {
      DocListResultVector results;
      idx->Search(tt_pattern, std::ref(results));
      results();
      return results;
    };

    return std::make_tuple(query, idx.size, idx.idx->GetSizeReport());
  };

  RegisterAllQueryBenchmarks(factory_result_vector_int,
                             factory.SequenceSize(),
                             idx_configs,
                             patterns_encoded,
                             query_bm_config,
                             FLAGS_min_s,
                             FLAGS_max_s,
                             2);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  //  benchmark::Shutdown();

  return 0;
}
