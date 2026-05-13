########################################################################
# CMake build script for DRet Library.

cmake_minimum_required(VERSION 3.10)

option(dret_build_tests "Build all of drl's tests." ON)
option(dret_build_tools "Build all of drl's tools." ON)
option(dret_build_benchmarks "Build all of drl's benchmarks." ON)

set(DRET_DEFAULT_ALPHABET_WIDTH 8 CACHE STRING "Default alphabet width for the ${PROJECT_NAME} project.")


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
set(CMAKE_CXX_STANDARD 20)


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
        ${PROJECT_SOURCE_DIR}/include/dret/doc_index_rmq_common.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_rmq.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_index_rmq.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_rmq_scheme.h
        ${PROJECT_SOURCE_DIR}/include/dret/rmq_get_doc_policies.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_gcda.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_freq_index_sada.h
        ${PROJECT_SOURCE_DIR}/include/dret/tf.h
        ${PROJECT_SOURCE_DIR}/include/dret/slp_tools.h
        ${PROJECT_SOURCE_DIR}/include/dret/algorithm.h
        ${PROJECT_SOURCE_DIR}/include/dret/run_length_encode.h
        ${PROJECT_SOURCE_DIR}/include/dret/config.h
        ${PROJECT_SOURCE_DIR}/include/dret/index_base.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_index.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_index_brute.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_sampled_tree.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_sampled_tree_gcda.h
        ${PROJECT_SOURCE_DIR}/include/dret/basic_slp_span_length.h
        ${PROJECT_SOURCE_DIR}/include/dret/differential_slp.h
        ${PROJECT_SOURCE_DIR}/include/dret/differential_light_slp.h
        ${PROJECT_SOURCE_DIR}/include/dret/doc_list_sampled_tree_dgcda.h

        # PDL (precomputed document listing) — append new include/dret/pdl/*.h here.
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/storage_policy.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/tree_core.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/tree_builder.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/build_pdl_core.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/set_codecs.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/get_docs.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/doc_list_pdl_plain.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/doc_list_pdl_rp.h
        ${PROJECT_SOURCE_DIR}/include/dret/pdl/doc_list_pdl_bc.h
)

add_library(dret INTERFACE)
target_sources(dret INTERFACE ${SOURCE_FILES})
target_include_directories(dret INTERFACE ${PROJECT_SOURCE_DIR}/include/dret)

add_compile_definitions(DRET_DEFAULT_ALPHABET_WIDTH=${DRET_DEFAULT_ALPHABET_WIDTH})

include(internal_utils)

include(ConfigSDSL)
include(ConfigGrammar)
include(ConfigSRIndex)
set(LIBS ${SDSL_LIB} grammar sr-index nlohmann_json::nlohmann_json)
target_link_libraries(dret INTERFACE ${LIBS})

include(ConfigRePair)
set(REPAIR_EXE "\"${${ExternalProjectName}_SOURCE_DIR}/bal/irepair\"")
target_compile_definitions(dret INTERFACE REPAIR_EXE=${REPAIR_EXE})

include(ConfigDSExtract)
set(VNMEXTRACT_EXE "\"${${ExternalProjectName}_VNMEXTRACT}\"")
target_compile_definitions(dret INTERFACE VNMEXTRACT_EXE=${VNMEXTRACT_EXE})


# TODO use find_package(gflags) instead find_library(...)
find_library(GFLAGS_LIB gflags)
if (NOT GFLAGS_LIB)
    include(ConfigGFlags)
    set(GFLAGS_LIB gflags::gflags)
endif ()

if (dret_build_tests)
    enable_testing()

    list(APPEND test_LIBS dret ${LIBS})

    include(ConfigGTest)
    list(APPEND test_LIBS gtest gtest_main gmock)

    cxx_test_with_flags_and_args(doc_freq_index_brute_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/doc_freq_index_brute_test.cpp)
    cxx_test_with_flags_and_args(tf_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/tf_test.cpp)
    cxx_test_with_flags_and_args(algorithm_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/algorithm_test.cpp)
    cxx_test_with_flags_and_args(doc_list_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/doc_list_test.cpp)
    cxx_test_with_flags_and_args(pdl_tree_builder_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/pdl_tree_builder_test.cpp)
    cxx_test_with_flags_and_args(pdl_tree_core_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/pdl_tree_core_test.cpp)
    cxx_test_with_flags_and_args(pdl_set_codecs_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/pdl_set_codecs_test.cpp)
    add_dependencies(pdl_set_codecs_test dsextract)
    cxx_test_with_flags_and_args(pdl_template_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/pdl_template_test.cpp)
    add_dependencies(pdl_template_test dsextract)
    cxx_test_with_flags_and_args(pdl_get_docs_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/pdl_get_docs_test.cpp)
    cxx_test_with_flags_and_args(pdl_index_test "" "${test_LIBS}" "" ${PROJECT_SOURCE_DIR}/test/pdl_index_test.cpp)
    add_dependencies(pdl_index_test dsextract)
endif ()


find_package(Threads)

if (dret_build_tools)
    include(ConfigRIndexLib)

    list(APPEND tool_LIBS dret ${LIBS} ${GFLAGS_LIB})

    find_package(Threads)
    list(APPEND tool_LIBS ${CMAKE_THREAD_LIBS_INIT})

    cxx_executable_with_flags(build_items "" "${tool_LIBS}" tool/build_items.cpp)
    cxx_executable_with_flags(build_items_sada "" "${tool_LIBS}" tool/build_items_sada.cpp tool/definitions.h)
    cxx_executable_with_flags(build_items_ilcp "" "${tool_LIBS}" tool/build_items_ilcp.cpp tool/definitions.h)
    cxx_executable_with_flags(build_items_gcda "" "${tool_LIBS}" tool/build_items_gcda.cpp tool/definitions.h)
    #    cxx_executable_with_flags(build_dslp "" "${tool_LIBS}" tool/build_dslp.cpp)
endif ()


if (dret_build_benchmarks)
    add_subdirectory(benchmark)
endif ()
