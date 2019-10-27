option(dret_build_tests "Build all of drl's tests." ON)

option(dret_build_tools "Build all of drl's tools." ON)

option(dret_build_benchmarks "Build all of drl's benchmarks." ON)


cmake_minimum_required(VERSION 2.8)


#Project Information
project(dret CXX)


# Location of internal cmake modules
set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${PROJECT_SOURCE_DIR}/cmake)

# Guard against in-source builds and bad build-type strings
include(ConfigSafeGuards)


if (dret_build_tests)
    include(ConfigGTest)
endif ()


if (dret_build_tools)
endif ()


if (dret_build_benchmarks)
endif ()

