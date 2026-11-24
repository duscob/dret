set(ExternalProjectName rlcsa)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/rlcsa.git
        GIT_TAG master
        FIND_PACKAGE_ARGS
)

set(RLCSA_ENABLE_BENCHMARKS OFF CACHE BOOL "")

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR})
