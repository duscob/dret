//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/10/26.
//
// PDL stored-set selection policy. Tiny enum-only header so both the
// query-time PDLTreeCore and the construction-time tree_builder can
// share it without one pulling in the other.
//

#pragma once

#include <cstdint>

namespace dret::pdl {

enum class StoragePolicy : uint8_t {
  // Default: a node is stored if it is childless, has the all-doc sentinel,
  // OR its pre-dedup occurrence weight exceeds storing_factor * |distinct
  // docs|. Ports drl's selection rule at drl/src/pdltree.cpp:390.
  OccurrenceWeighted = 0,
  // Diagnostic upper bound: every node stored, regardless of weight.
  StoreAllInternal = 1,
  // Ablation: only childless (collapsed-block + explicit-leaf) nodes
  // stored; internal nodes kept in the tree for navigation but not
  // contributing a stored set.
  LeavesOnly = 2,
};

}  // namespace dret::pdl
