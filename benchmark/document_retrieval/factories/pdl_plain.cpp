//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the PDL cells with plain stored sets, compiled on their own (see cells.h).
//

#include "pdl.h"

namespace bench::factories::pdl {

MakeIndex MakePlain(const Params& t_params) {
  return ByGetDoc<IdxPlain>(t_params);
}

}  // namespace bench::factories::pdl
