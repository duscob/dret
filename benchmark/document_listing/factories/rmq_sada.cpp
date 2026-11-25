//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the RMQ cells with the SADA (Sadakane, marking) core, compiled on their own (see cells.h).
//

#include "rmq.h"

namespace bench::factories::rmq {

namespace {
template <typename TGetDoc>
using Core = SadaCore<Storage, TGetDoc>;
}  // namespace

MakeIndex MakeSada(const Params& t_params) {
  return ByGetDoc<Core>(t_params);
}

}  // namespace bench::factories::rmq
