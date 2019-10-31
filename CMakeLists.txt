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
set(CMAKE_CXX_STANDARD 17)


include_directories(${PROJECT_SOURCE_DIR}/include)

set(SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/include/dret/tf.h)

add_library(dret INTERFACE)
target_sources(dret INTERFACE ${SOURCE_FILES})
target_include_directories(dret INTERFACE ${PROJECT_SOURCE_DIR}/include/dret)

include(internal_utils)

if (dret_build_tests)
    enable_testing()

    include(ConfigGTest)

    find_library(SDSL_LIB sdsl REQUIRED)

    cxx_test_with_flags_and_args(tf_test "" "gtest;gtest_main;gmock;sdsl" "" ${PROJECT_SOURCE_DIR}/test/tf_test.cpp)
endif ()


if (dret_build_tools)
endif ()


if (dret_build_benchmarks)
endif ()

