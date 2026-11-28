//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the RMQ cells with the ILCP-V (value rule) core, compiled on their own (see cells.h).
//

#include "rmq.h"

namespace bench::factories::rmq {

namespace {
// The core over a backend, its stored values (run_values) in TValues.
template <typename TValues>
struct WithValues {
  template <typename TGetDoc>
  using Core = IlcpValueCore<Storage, TGetDoc, TValues>;
};
}  // namespace

MakeIndex MakeIlcpValue(const Params& t_params) {
  switch (t_params.run_values) {
    case axes::RunValuesVariant::IV:
      return ByGetDoc<WithValues<RunValues_IV>::Core>(t_params);
    case axes::RunValuesVariant::VV:
      return ByGetDoc<WithValues<RunValues_VV>::Core>(t_params);
    case axes::RunValuesVariant::DV:
    default:
      return ByGetDoc<WithValues<RunValues_DV>::Core>(t_params);
  }
}

}  // namespace bench::factories::rmq
