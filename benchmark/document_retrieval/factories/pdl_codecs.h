//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The PDL cells' index makers, one per stored-set codec, each defined in a source file of its own
// (pdl_<codec>.cpp; see cells.h for why). pdl.cpp turns the sweep into cells and calls these.
//

#pragma once

#include <cstdint>

#include "../axes.h"
#include "cells.h"

namespace bench::factories::pdl {

// One PDL cell.
struct Params {
  std::uint32_t block_size;
  float storing_factor;
  axes::PDLStoragePolicy policy;
  axes::GetDocEnum get_doc;
  axes::GCDASLPVariant gcda_slp;
  axes::BareSLPVariant bare_slp;
  axes::DGCDASLPVariant dgcda_slp;
  // Block size of a backend built apart from PDL's own (GCDA-light, a differential bare grammar).
  std::uint32_t backend_block_size;
};

MakeIndex MakePlain(const Params& t_params);  // pdl_plain.cpp
MakeIndex MakeRP(const Params& t_params);     // pdl_rp.cpp
MakeIndex MakeBC(const Params& t_params);     // pdl_bc.cpp

}  // namespace bench::factories::pdl
