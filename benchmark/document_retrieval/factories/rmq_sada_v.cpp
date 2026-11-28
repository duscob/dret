//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the RMQ cells with the SADA-V (value rule) core, compiled on their own (see cells.h).
//

#include "rmq.h"

namespace bench::factories::rmq {

namespace {
// The core over a backend, its stored values (prev_doc) in TValues.
template <typename TValues>
struct WithValues {
  template <typename TGetDoc>
  using Core = SadaValueCore<Storage, TGetDoc, TValues>;
};
}  // namespace

MakeIndex MakeSadaValue(const Params& t_params) {
  switch (t_params.prev_doc) {
    case axes::PrevDocVariant::DV:
      return ByGetDoc<WithValues<PrevDoc_DV>::Core>(t_params);
    case axes::PrevDocVariant::VV:
      return ByGetDoc<WithValues<PrevDoc_VV>::Core>(t_params);
    case axes::PrevDocVariant::IV:
    default:
      return ByGetDoc<WithValues<PrevDoc_IV>::Core>(t_params);
  }
}

}  // namespace bench::factories::rmq
