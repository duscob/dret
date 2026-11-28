//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode factory of GCDA over a differential grammar, compiled on its own (see query.h).
//

#include "dgcda.h"
#include "query.h"

namespace bench::factories::query {

Built DGCDA(Storage t_storage,
            dret::Config& t_config,
            uint32_t t_block_size,
            float t_storing_factor,
            axes::DGCDASLPVariant t_slp,
            axes::GCDASetsCodec t_sets) {
  return dgcda::Make(t_storage, t_config, t_block_size, t_storing_factor, t_slp, t_sets);
}

}  // namespace bench::factories::query
