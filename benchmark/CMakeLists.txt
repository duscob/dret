# Configure benchmarks
set(benchmark_LIBS dret ${LIBS} gflags::gflags)

if (NOT TARGET benchmark)
    include(ConfigGBenchmark)
endif ()
list(APPEND benchmark_LIBS benchmark::benchmark benchmark::benchmark_main)

find_package(Threads)
list(APPEND benchmark_LIBS ${CMAKE_THREAD_LIBS_INIT})

#    cxx_executable_with_flags(query_doc_freq_index_bm "" "${benchmark_LIBS}" query_doc_freq_index_bm.cpp)
cxx_executable_with_flags(bm_query_doc_freq_index "" "${benchmark_LIBS}" document_frequencies/bm_query_doc_freq_index.cpp)

cxx_executable_with_flags(bm_build_items "" "${benchmark_LIBS}" document_listing/bm_build_items.cpp)
cxx_executable_with_flags(bm_query_doc_list "" "${benchmark_LIBS}" document_listing/bm_query_doc_list.cpp)
cxx_executable_with_flags(bm_doc_list "" "${benchmark_LIBS}" document_listing/bm_doc_list.cpp)
