//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 5/14/26.
//
// Header inclusion / two-phase lookup canary, Task 46 of
// docs/pdl_indexes_tasks.md.
//
// Includes ONLY the public PDL Plain header — no extra transitive
// helpers — and forces compile-time resolution of the call sites
// that depend on declarations being visible at template-definition
// point:
//
//   - PDLTreeCore<>::computeCoverFull        (Task 12 protected virtual,
//                                              called unqualified from
//                                              DLSampledTreeScheme::Search)
//   - dret::pdl::collectSizes(SizeReport&,
//                              PlainCodec<>&,
//                              prefix)         (Task 32 free-function
//                                              overload, ADL on the codec
//                                              type)
//   - dret::pdl::collectSizes(SizeReport&,
//                              PDLTreeCore<>&,
//                              prefix)         (Task 32 free-function
//                                              overload, ADL on the core
//                                              type)
//
// If a future refactor moves any of these declarations OUT of the
// transitive include set of doc_list_pdl_plain.h (e.g., into a
// separate header not pulled in by the public PDL Plain header), or
// declares them AFTER the class-template definition that calls them
// unqualified, this test fails to compile. The runtime assertions
// are trivial — the load-bearing assertion is "this file compiles".
//

// Deliberately the only PDL include — no other transitive helpers.
#include "dret/doc_list/doc_list_pdl_plain.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace {

// 1. The class template itself must instantiate at default template
//    arguments. This pulls in the full member-function set, including
//    overrides of computeCoverFull / GetSizeReport / loadInner, and
//    forces the compiler to parse them. If any refactor breaks
//    two-phase lookup of the helpers those overrides depend on, the
//    instantiation fails here.
TEST(PDLHeaderInclusion, DocListIdxPDLPlainInstantiates) {
  static_assert(sizeof(dret::pdl::DocListIdxPDLPlain<>) > 0);
  dret::pdl::DocListIdxPDLPlain<> idx;
  (void)idx;
}

// 2. PDLTreeCore::computeCoverFull is what DLSampledTreeScheme's
//    Search() forwards to via PDL Plain's override. Force its
//    instantiation at default template args; an empty core returns
//    no raw ranges and no nodes.
TEST(PDLHeaderInclusion, PDLTreeCoreComputeCoverFullCompilesAndRuns) {
  dret::pdl::PDLTreeCore<> core;
  std::vector<std::pair<std::size_t, std::size_t>> raw_ranges;
  std::vector<std::size_t> nodes;
  core.computeCoverFull(/*sp=*/0, /*ep=*/0, raw_ranges, nodes);
  EXPECT_THAT(raw_ranges, testing::IsEmpty());
  EXPECT_THAT(nodes, testing::IsEmpty());
}

// 3. ADL on PlainCodec<> must find dret::pdl::collectSizes. Calling
//    unqualified is the load-bearing form: ADL is what surfaces the
//    overload from the codec's namespace at template-definition
//    point. If set_codecs.h's collectSizes overload is moved out of
//    PDL Plain's transitive include set, the call fails to compile.
TEST(PDLHeaderInclusion, CollectSizesPlainCodecResolvesViaADL) {
  dret::pdl::PlainCodec<> codec;
  dret::SizeReport report;
  collectSizes(report, codec, "p_");  // ADL on dret::pdl::PlainCodec<>
  // PlainCodec::GetSizeReport emits one row per stored container
  // even when empty (chunks_objs, chunks_pos); we don't pin specific
  // bytes — the compile-time resolution is the assertion. Verify
  // the prefix at least made it onto every reported field.
  for (const auto& f : report) {
    EXPECT_NE(f.name.find("p_"), std::string::npos) << f.name;
  }
}

// 4. Same canary for PDLTreeCore: unqualified collectSizes must
//    resolve via ADL on dret::pdl::PDLTreeCore<>. tree_core.h's
//    overload must remain in scope.
TEST(PDLHeaderInclusion, CollectSizesPDLTreeCoreResolvesViaADL) {
  dret::pdl::PDLTreeCore<> core;
  dret::SizeReport report;
  collectSizes(report, core, "c_");  // ADL on dret::pdl::PDLTreeCore<>
  // An un-Assembled core's GetSizeReport returns rows for every
  // empty bitvector / int_vector member; we don't pin specific
  // bytes — the compile-time resolution is the assertion.
  (void)report;
}

}  // namespace
