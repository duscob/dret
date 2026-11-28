//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the RMQ cells with the ILCP (Gagie et al., marking) core, compiled on their own (see cells.h).
//

#include "rmq.h"

namespace bench::factories::rmq {

namespace {
template <typename TGetDoc>
using Core = IlcpCore<Storage, TGetDoc>;
}  // namespace

MakeIndex MakeIlcp(const Params& t_params) {
  return ByGetDoc<Core>(t_params);
}

}  // namespace bench::factories::rmq
