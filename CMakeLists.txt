option(dret_build_tests "Build all of drl's tests." ON)

option(dret_build_tools "Build all of drl's tools." ON)

option(dret_build_benchmarks "Build all of drl's benchmarks." ON)


cmake_minimum_required(VERSION 2.8)


#Project Information
project(dret CXX)


if (dret_build_tests)
endif ()


if (dret_build_tools)
endif ()


if (dret_build_benchmarks)
endif ()

