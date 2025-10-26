set(ExternalProjectName r-index)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/r-index.git
        GIT_TAG library
        FIND_PACKAGE_ARGS
)

set(rindex_build_tools ON CACHE BOOL "")
set(rindex_build_tests OFF CACHE BOOL "")
set(rindex_build_benchmarks OFF CACHE BOOL "")
set(rindex_install OFF CACHE BOOL "")

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/internal)
