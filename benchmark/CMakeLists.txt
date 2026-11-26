# Configure benchmarks
set(benchmark_LIBS dret ${LIBS} gflags::gflags)

if (NOT TARGET benchmark)
    include(ConfigGBenchmark)
endif ()
list(APPEND benchmark_LIBS benchmark::benchmark benchmark::benchmark_main)

find_package(Threads)
list(APPEND benchmark_LIBS ${CMAKE_THREAD_LIBS_INIT})

cxx_executable_with_flags(bm_query_doc_freq_index "" "${benchmark_LIBS}" document_frequencies/bm_query_doc_freq_index.cpp)

# bm_doc_list: each family's cells (document_listing/factories/cells.h) are
# compiled in a source file of its own, RMQ one per core and PDL one per codec,
# so a change to one family does not move another's query times through GCC's
# inlining choices (cells.h says how it did). bm_doc_list.cpp is the driver.
set(bm_doc_list_family_sources
        document_listing/factories/brute.cpp
        document_listing/factories/gcda.cpp
        document_listing/factories/gcda_backend.cpp
        document_listing/factories/dgcda.cpp
        document_listing/factories/slp_ns.cpp
        document_listing/factories/rmq.cpp
        document_listing/factories/rmq_sada.cpp
        document_listing/factories/rmq_ilcp.cpp
        document_listing/factories/rmq_cilcp.cpp
        document_listing/factories/rmq_sada_v.cpp
        document_listing/factories/rmq_ilcp_v.cpp
        document_listing/factories/rmq_cilcp_v.cpp
        document_listing/factories/pdl.cpp
        document_listing/factories/pdl_plain.cpp
        document_listing/factories/pdl_rp.cpp
        document_listing/factories/pdl_bc.cpp)
cxx_executable_with_flags(bm_doc_list "" "${benchmark_LIBS}" ${bm_doc_list_family_sources} document_listing/bm_doc_list.cpp)
cxx_executable_with_flags(migrate_cache_components "" "${benchmark_LIBS}" document_listing/migrate_cache_components.cpp)

# The benchmarks shell out to two external-project binaries at run time: PDL-BC
# to vnmextract (dsextract) and GCDA/DGCDA to irepair (repair). Their paths are
# baked in as VNMEXTRACT_EXE / REPAIR_EXE_BAL32 / REPAIR_EXE_BAL64 compile
# definitions on the `dret` INTERFACE target, but a compile definition creates
# no build ordering, and an INTERFACE library cannot carry the dependency for
# its consumers.
#
# Without this, `cmake --build <dir> --target bm_doc_list` in a fresh tree
# yields a binary whose VNMEXTRACT_EXE points at a file that was never built.
# PDL-BC then reports zero bicliques on every collection -- silently, and
# indistinguishably from a genuine "nothing to compress" result. That cost a
# long misdiagnosis in 2026-08. The test targets already declared this
# dependency, which is exactly why the tests passed while the benchmarks were
# quietly broken.
foreach (_dret_bm bm_query_doc_freq_index bm_doc_list)
    if (TARGET ${_dret_bm})
        add_dependencies(${_dret_bm} dsextract repair)
    endif ()
endforeach ()
