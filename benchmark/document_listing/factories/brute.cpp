//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode factory of the brute-force baselines and CountRI, compiled on its own (see query.h).
//

#include "brute.h"
#include "query.h"

namespace bench::factories::query {

// The r-index and sr-index baselines only load: their caches come from a construct run, and
// DocListIdxBrute's construct() would rebuild them from the text.
Built BruteR(Storage t_storage, dret::Config& t_config) {
  auto idx = std::make_shared<brute::IdxR<Storage>>(t_storage);
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

Built BruteSr(Storage t_storage, dret::Config& t_config, std::size_t t_sampling) {
  auto idx = std::make_shared<brute::IdxSr<Storage>>(t_storage, sri::SrIndexValidArea<Storage>(t_storage, t_sampling));
  idx->load(t_config);
  return {idx, sdsl::size_in_bytes(*idx)};
}

Built BruteDA(Storage t_storage, dret::Config& t_config) {
  return brute::MakeBackend<brute::IdxDA<Storage>>(t_storage, t_config);
}

Built BruteRLCSA(Storage t_storage, dret::Config& t_config) {
  return brute::MakeBackend<brute::IdxRLCSA<Storage>>(t_storage, t_config);
}

Built CountRI(Storage t_storage, dret::Config& t_config) {
  return brute::MakeBackend<brute::IdxCount<Storage>>(t_storage, t_config);
}

// Loads the r-index on its first call, so it needs the same cache files as Brute-RI.
std::function<std::size_t(const std::string&)> OccurrenceCounter(Storage t_storage, dret::Config& t_config) {
  std::shared_ptr<sri::RIndex<Storage>> r_index;
  return [t_storage, t_config, r_index](const std::string& t_pattern) mutable -> std::size_t {
    if (!r_index) {
      r_index = std::make_shared<sri::RIndex<Storage>>(t_storage);
      r_index->load(t_config);
    }
    const auto [start, end] = r_index->Count(t_pattern);
    return end > start ? end - start : 0;
  };
}

}  // namespace bench::factories::query
