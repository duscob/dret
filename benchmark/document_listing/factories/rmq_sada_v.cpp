//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// Query-mode factory of the RMQ index with the SADA-V (value rule) core, compiled on its own (see query.h).
//

#include "query.h"
#include "rmq.h"

namespace bench::factories::query {

Built SadaValue(Storage t_storage, dret::Config& t_config, const RMQArgs& t_args) {
  return rmq::detail::MakeOneSadaValue<rmq::SadaValueCore>(t_storage,
                                                           t_config,
                                                           t_args.block_size,
                                                           t_args.storing_factor,
                                                           t_args.get_doc,
                                                           t_args.gcda_slp,
                                                           t_args.bare_slp,
                                                           t_args.dgcda_slp,
                                                           t_args.prev_doc,
                                                           t_args.sa_sampling);
}

}  // namespace bench::factories::query
