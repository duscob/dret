//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode factory of PDL with plain stored sets, compiled on its own (see query.h).
//

#include "pdl.h"
#include "query.h"

namespace bench::factories::query {

Built PDLPlain(Storage t_storage, dret::Config& t_config, const PDLArgs& t_args) {
  return pdl::detail::MakeForCodec<pdl::IdxPlain>(t_storage,
                                                  t_config,
                                                  t_args.block_size,
                                                  t_args.storing_factor,
                                                  axes::toPDLStoragePolicy(t_args.policy),
                                                  t_args.get_doc,
                                                  t_args.gcda_slp,
                                                  t_args.bare_slp,
                                                  t_args.dgcda_slp,
                                                  t_args.sa_sampling,
                                                  t_args.backend_block_size);
}

}  // namespace bench::factories::query
