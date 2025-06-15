########################################################################
# CMake build script for DRet Library.

cmake_minimum_required(VERSION 3.10)

option(dret_build_tests "Build all of drl's tests." ON)

option(dret_build_tools "Build all of drl's tools." ON)

option(dret_build_benchmarks "Build all of drl's benchmarks." ON)


#Project Information
project(dret
        VERSION 0.0.0
        LANGUAGES CXX)


# Location of internal cmake modules
set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${PROJECT_SOURCE_DIR}/cmake)

# Guard against in-source builds and bad build-type strings
include(ConfigSafeGuards)


#Global Setup
set(CMAKE_CXX_STANDARD 17)


#Global C++ Flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")


# Set common include folder for module
set(COMMON_INCLUDES
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_INSTALL_PREFIX}/include
        ${CMAKE_PREFIX_PATH}/include)

include_directories(${COMMON_INCLUDES})


set(SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_brute.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_rmq.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_gcda.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_sada.h
        ${PROJECT_SOURCE_DIR}/include/dret/tf.h
        ${PROJECT_SOURCE_DIR}/include/dret/algorithm.h
        ${PROJECT_SOURCE_DIR}/include/dret/run_length_encode.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_index.h)

add_library(dret INTERFACE)
target_sources(dret INTERFACE ${SOURCE_FILES})
target_include_directories(dret INTERFACE ${PROJECT_SOURCE_DIR}/include/dret)

include(internal_utils)

find_package(sdsl 3.0)
if (NOT sdsl_FOUND)
    include(ConfigSDSL)
endif ()

# TODO use find_package(gflags) instead find_library(...)
find_library(GFLAGS_LIB gflags)
if (NOT GFLAGS_LIB)
    include(ConfigGFlags)
    set(GFLAGS_LIB gflags::gflags)
endif ()


#find_library(GRAMMAR_LIB grammar REQUIRED)
include(ConfigGrammar)
set(GRAMMAR_LIB grammar)
set(LIBS ${SDSL_LIB} ${GRAMMAR_LIB})

if (dret_build_tests)
    enable_testing()

    include(ConfigGTest)

    cxx_test_with_flags_and_args(doc_freq_index_brute_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/doc_freq_index_brute_test.cpp)
    cxx_test_with_flags_and_args(tf_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/tf_test.cpp)
    cxx_test_with_flags_and_args(algorithm_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/algorithm_test.cpp)
    cxx_test_with_flags_and_args(doc_list_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/doc_list_test.cpp)
endif ()

find_package(Threads)

if (dret_build_tools)
    include(ConfigRIndexLib)
    include(ConfigNJSON)
    include(ConfigSRIndex)

    cxx_executable_with_flags(build_items "" "${SDSL_LIBS};${GFLAGS_LIB};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items.cpp)
    cxx_executable_with_flags(build_items_sada "" "${SDSL_LIBS};${GFLAGS_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items_sada.cpp tool/definitions.h)
    cxx_executable_with_flags(build_items_ilcp "" "${SDSL_LIBS};${GFLAGS_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items_ilcp.cpp tool/definitions.h)
    cxx_executable_with_flags(build_items_gcda "" "${SDSL_LIBS};${GFLAGS_LIB};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items_gcda.cpp tool/definitions.h)
    #    cxx_executable_with_flags(build_dslp "" "${SDSL_LIBS};${GFLAGS_LIB};${GRAMMAR_LIB}" tool/build_dslp.cpp)
endif ()


if (dret_build_benchmarks)
    if (NOT TARGET benchmark)
        include(ConfigGBenchmark)
    endif ()

    #    cxx_executable_with_flags(query_doc_freq_index_bm "" "${GFLAGS_LIB};benchmark;dret;${LIBS};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" benchmark/query_doc_freq_index_bm.cpp)
    cxx_executable_with_flags(bm_query_doc_freq_index "" "${GFLAGS_LIB};benchmark;dret;${LIBS};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" benchmark/bm_query_doc_freq_index.cpp)

    cxx_executable_with_flags(bm_query_doc_list "" "${GFLAGS_LIB};benchmark;dret;${LIBS};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" benchmark/document_listing/bm_query_doc_list.cpp)
endif ()

