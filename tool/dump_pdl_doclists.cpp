//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/14/26.
//
// Standalone harness for Task 44 (drl parity smoke). Builds dret
// PDL Plain (DA backing) over a collection and dumps per-pattern
// document sets to stdout in the same format as drl's
// `build_docs_list` (RLCSA brute-force):
//
//   <pattern>\n
//   <doc> <doc> ... \n
//
// Doc ids are zero-based, sorted ascending. The output of this tool
// against any small dataset can be diffed directly against drl's
// build_docs_list output (RLCSA brute-force, also zero-based on the
// drl smoke fixture) to cross-check that dret PDL agrees with the
// independently-implemented drl reference.
//
// drl's PDLRP is in turn validated against drl's own RLCSA brute
// force in the drl test suite, so agreement here transitively
// validates dret PDL against drl PDLRP.
//

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <gflags/gflags.h>

#include "dret/config.h"
#include "dret/doc_list/doc_list_pdl.h"

DEFINE_string(data, "", "Data file (collection). Document delimiter is byte 0x00.");
DEFINE_string(patterns, "", "Patterns file (one pattern per line).");
DEFINE_string(data_dir, "", "Cache directory. Defaults to dirname(--data).");
DEFINE_int32(block_size, 512, "PDL block size.");
DEFINE_double(storing_factor, 4.0, "PDL storing factor.");

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_data.empty() || FLAGS_patterns.empty()) {
    std::cerr << "Usage: dump_pdl_doclists --data=<file> --patterns=<file>\n";
    return 1;
  }

  const std::filesystem::path data_path = FLAGS_data;
  const std::filesystem::path data_dir = FLAGS_data_dir.empty()
      ? data_path.parent_path()
      : std::filesystem::path(FLAGS_data_dir);

  // Build PDL Plain via the standard construct() pipeline.
  dret::Config config(data_path, data_dir, sri::SDSL_LIBDIVSUFSORT,
                      /*delete_files=*/false, /*data_width=*/8, /*data_delim=*/0);

  dret::GenericStorage storage;
  dret::pdl::DocListIdxPDLPlain<dret::GenericStorage> pdl(
      storage,
      static_cast<uint32_t>(FLAGS_block_size),
      static_cast<float>(FLAGS_storing_factor));
  construct(pdl, config);

  // Read patterns line-by-line, query, dump in drl's build_docs_list
  // format (pattern, doc-list, blank line for the final-newline padding
  // drl emits).
  std::ifstream pattern_file(FLAGS_patterns);
  if (!pattern_file) {
    std::cerr << "Cannot open patterns file: " << FLAGS_patterns << "\n";
    return 2;
  }

  std::string line;
  while (std::getline(pattern_file, line)) {
    if (line.empty()) continue;
    std::vector<std::size_t> docs;
    auto report = [&docs](std::size_t d) { docs.push_back(d); };
    pdl.Search(line, report);
    std::sort(docs.begin(), docs.end());
    docs.erase(std::unique(docs.begin(), docs.end()), docs.end());

    std::cout << line << "\n";
    for (auto d : docs) std::cout << d << " ";
    std::cout << "\n";
  }
  std::cout << "\n";  // drl emits one trailing blank line; match it.
  return 0;
}
