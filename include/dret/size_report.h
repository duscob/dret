//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/24/26.
//

#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/io.hpp>

#include <grammar/sampled_slp.h>

namespace dret {

struct SizeField {
  std::string name;
  std::size_t bytes = 0;
};

using SizeReport = std::vector<SizeField>;

inline void append(SizeReport& report, std::string name, std::size_t bytes) {
  report.push_back({std::move(name), bytes});
}

inline std::size_t totalBytes(const SizeReport& report) {
  std::size_t total = 0;
  for (const auto& f : report)
    total += f.bytes;
  return total;
}

inline void writeSizesJson(const std::string& path, const SizeReport& report) {
  std::ofstream out(path);
  out << "{\n";
  out << "  \"total_bytes\": " << totalBytes(report) << ",\n";
  out << "  \"fields\": {\n";
  for (std::size_t i = 0; i < report.size(); ++i) {
    out << "    \"" << report[i].name << "\": " << report[i].bytes;
    if (i + 1 < report.size())
      out << ",";
    out << "\n";
  }
  out << "  }\n";
  out << "}\n";
}

//~~~~~~~ collectSizes overloads ~~~~~~~

// grammar::LightSLP<_SLP, _SampledSLP, _Chunks> (GCDA inner type)
template <typename TSLP, typename TSampledSLP, typename TChunks>
void collectSizes(SizeReport& out,
                  const grammar::LightSLP<TSLP, TSampledSLP, TChunks>& slp,
                  const std::string& prefix = "") {
  append(out, prefix + "base_slp", sdsl::size_in_bytes(static_cast<const TSLP&>(slp)));
  append(out, prefix + "sampled_slp", sdsl::size_in_bytes(static_cast<const TSampledSLP&>(slp)));
  append(out, prefix + "covers", sdsl::size_in_bytes(slp.GetCovers()));
}

}  // namespace dret
