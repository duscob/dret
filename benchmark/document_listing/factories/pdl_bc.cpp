//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The indexes of the PDL cells with biclique-compressed stored sets, compiled on their own (see cells.h).
//

#include "pdl.h"

namespace bench::factories::pdl {

MakeIndex MakeBC(const Params& t_params) {
  return ByGetDoc<IdxBC>(t_params);
}

}  // namespace bench::factories::pdl
