option(dret_build_tests "Build all of drl's tests." ON)

option(dret_build_tools "Build all of drl's tools." ON)

option(dret_build_benchmarks "Build all of drl's benchmarks." ON)


cmake_minimum_required(VERSION 3.0)


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
set(CMAKE_CXX_STANDARD 14)


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
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_sada.h
        ${PROJECT_SOURCE_DIR}/include/dret/tf.h
        ${PROJECT_SOURCE_DIR}/include/dret/algorithm.h)

add_library(dret INTERFACE)
target_sources(dret INTERFACE ${SOURCE_FILES})
target_include_directories(dret INTERFACE ${PROJECT_SOURCE_DIR}/include/dret)

include(internal_utils)

find_library(SDSL_LIB sdsl)
find_library(DIVSUFSORT_LIB divsufsort)
find_library(DIVSUFSORT64_LIB divsufsort64)
set(SDSL_LIBS ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB})
find_library(GRAMMAR_LIB grammar REQUIRED)
set(LIBS ${SDSL_LIB} ${GRAMMAR_LIB})

if (dret_build_tests)
    enable_testing()

    include(ConfigGTest)

    cxx_test_with_flags_and_args(doc_freq_index_brute_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/doc_freq_index_brute_test.cpp)
    cxx_test_with_flags_and_args(tf_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/tf_test.cpp)
    cxx_test_with_flags_and_args(algorithm_test "" "gtest;gtest_main;gmock;dret;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/algorithm_test.cpp)
endif ()

# TODO use find_package(gflags) instead find_library(...)
find_library(GFLAGS_LIB gflags)
find_package(Boost COMPONENTS filesystem)
find_package(Threads)

if (dret_build_tools)
    cxx_executable_with_flags(build_items "" "${SDSL_LIBS};${GFLAGS_LIB};${Boost_LIBRARIES};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items.cpp)
    cxx_executable_with_flags(build_items_sada "" "${SDSL_LIBS};${GFLAGS_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items_sada.cpp tool/definitions.h)
    cxx_executable_with_flags(build_items_ilcp "" "${SDSL_LIBS};${GFLAGS_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items_ilcp.cpp tool/definitions.h)
    cxx_executable_with_flags(build_items_gcda "" "${SDSL_LIBS};${GFLAGS_LIB};${Boost_LIBRARIES};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" tool/build_items_gcda.cpp tool/definitions.h)
#    cxx_executable_with_flags(build_dslp "" "${SDSL_LIBS};${GFLAGS_LIB};${GRAMMAR_LIB}" tool/build_dslp.cpp)
endif ()


if (dret_build_benchmarks)
    include(ConfigGBenchmark)

#    cxx_executable_with_flags(query_doc_freq_index_bm "" "${GFLAGS_LIB};benchmark;dret;${LIBS};${Boost_LIBRARIES};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" benchmark/query_doc_freq_index_bm.cpp)
    cxx_executable_with_flags(bm_query_doc_freq_index "" "${GFLAGS_LIB};benchmark;dret;${LIBS};${Boost_LIBRARIES};${GRAMMAR_LIB};${CMAKE_THREAD_LIBS_INIT}" benchmark/bm_query_doc_freq_index.cpp)
endif ()

