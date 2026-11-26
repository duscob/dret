//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/25/26.
//
// The RMQ cells' index makers, one per core, each defined in a source file of its own
// (rmq_<core>.cpp; see cells.h for why). rmq.cpp turns the sweep into cells and calls these.
//

#pragma once

#include <cstdint>

#include "../axes.h"
#include "cells.h"

namespace bench::factories::rmq {

// One RMQ cell. The published cores read neither run_values nor prev_doc; SADA-V reads prev_doc,
// ILCP-V and CILCP-V run_values.
struct Params {
  std::uint32_t block_size;  // of a grammar backend (a differential bare one: its sample spacing)
  float storing_factor;      // of a sampled grammar backend
  axes::GetDocEnum get_doc;
  axes::GCDASLPVariant gcda_slp;
  axes::BareSLPVariant bare_slp;
  axes::DGCDASLPVariant dgcda_slp;
  axes::RunValuesVariant run_values;
  axes::PrevDocVariant prev_doc;
};

MakeIndex MakeSada(const Params& t_params);        // rmq_sada.cpp
MakeIndex MakeIlcp(const Params& t_params);        // rmq_ilcp.cpp
MakeIndex MakeCilcp(const Params& t_params);       // rmq_cilcp.cpp
MakeIndex MakeSadaValue(const Params& t_params);   // rmq_sada_v.cpp
MakeIndex MakeIlcpValue(const Params& t_params);   // rmq_ilcp_v.cpp
MakeIndex MakeCilcpValue(const Params& t_params);  // rmq_cilcp_v.cpp

}  // namespace bench::factories::rmq
