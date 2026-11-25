//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode factory of Brute-GCDA (brute force over a grammar without a tree), compiled on its own (see query.h).
//

#include "slp_ns.h"
#include "query.h"

namespace bench::factories::query {

Built SLPNS(Storage t_storage,
            dret::Config& t_config,
            axes::BareSLPVariant t_slp,
            uint32_t t_spacing,
            std::optional<axes::GCDASLPVariant> t_sampled,
            float t_storing_factor) {
  return slp_ns::Make(t_storage, t_config, t_slp, t_spacing, t_sampled, t_storing_factor);
}

}  // namespace bench::factories::query
