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
  // The rule of Gagie et al.'s PDL (drl/src/pdltree.cpp, with storeThisSet):
  // a node is stored if childless, has the all-doc sentinel, or the lists a
  // query would read instead -- those of its nearest stored descendants --
  // total more than storing_factor * |its distinct docs|. OccurrenceWeighted
  // counts occurrences instead of those lists, so it stores more nodes.
  ListWeighted = 3,
};

// Whether a policy's selection reads the storing factor. Only the
// occurrence-weighted rule does; a core built under another policy is the same
// index for every storing factor, and is cached and reported once.
constexpr bool PolicyReadsStoringFactor(StoragePolicy t_policy) {
  return t_policy == StoragePolicy::OccurrenceWeighted || t_policy == StoragePolicy::ListWeighted;
}

}  // namespace dret::pdl
