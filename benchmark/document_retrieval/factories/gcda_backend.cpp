//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The cells of GCDA over a document-array backend: the sampled tree and node lists of the (block
// size, storing factor) cell, the ends of a range expanded by the backend (DA, SA-Phi or RLCSA)
// instead of the grammar.
//

#include <stdexcept>
#include <string>
#include <vector>

#include "dret/backend/da.h"
#include "dret/backend/rlcsa.h"
#include "dret/backend/sa_phi.h"
#include "dret/doc_list/doc_list_gcda_backend.h"

#include "cells.h"
#include "gcda.h"
#include "names.h"

namespace bench::factories {

namespace {

using namespace bench::axes;

template <typename TGetDocs, typename TSets>
using Idx = dret::gcda::
    DocListIdxGCDABackend<Storage, dret::Alphabet<>, sri::RIndexCount<Storage, dret::Alphabet<>>, TGetDocs, TSets>;

template <typename TSets>
MakeIndex Make(GetDocEnum t_get_doc, std::uint32_t t_bs, float t_sf) {
  constexpr auto kW = dret::Alphabet<>::int_width;
  switch (t_get_doc) {
    case GetDocEnum::SAPhiR:
      return MakeOf<Idx<dret::backend::SAPhi<Storage, kW>, TSets>>(t_bs, t_sf);
    case GetDocEnum::RLCSA:
      return MakeOf<Idx<dret::backend::RLCSA<Storage, kW>, TSets>>(t_bs, t_sf);
    case GetDocEnum::DA:
      return MakeOf<Idx<dret::backend::PlainDA<Storage, kW>, TSets>>(t_bs, t_sf);
    default:
      throw std::invalid_argument("gcda_backend: get_doc must be da, sa_phi_r or rlcsa");
  }
}

}  // namespace

std::vector<Cell> Cells(const spec::GCDABackendSweep& t_sweep) {
  using namespace bench::names;
  std::vector<Cell> cells;
  for (auto gd : t_sweep.get_doc)
    for (auto sets : t_sweep.sets)
      for (auto bs : t_sweep.block_size)
        for (auto sf : t_sweep.storing_factor) {
          Cell cell{.name = KeyedName("GCDA",
                                      {{"tree", "sampled"},
                                       {"get-doc", GetDocValue(gd)},
                                       {"sets", SetsValue(sets)},
                                       {"block-size", IntStr(bs)},
                                       {"storing-factor", SfStr(sf)}}),
                    .make = sets == GCDASetsCodec::Plain ? Make<gcda::Sets_Plain>(gd, bs, sf)
                                                         : Make<gcda::Sets_RP>(gd, bs, sf)};
          cell.block_size = bs;
          cell.storing_factor = sf;
          cells.push_back(std::move(cell));
        }
  return cells;
}

}  // namespace bench::factories
