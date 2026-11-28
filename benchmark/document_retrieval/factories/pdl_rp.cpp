//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the PDL cells with Re-Pair stored sets, compiled on their own (see cells.h).
//

#include "pdl.h"

namespace bench::factories::pdl {

MakeIndex MakeRP(const Params& t_params) {
  return ByGetDoc<IdxRP>(t_params);
}

}  // namespace bench::factories::pdl
