set(ExternalProjectName grammar)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/grammar.git
        GIT_TAG compact
        FIND_PACKAGE_ARGS
)

set(${ExternalProjectName}_build_tools ON CACHE BOOL "")
set(${ExternalProjectName}_build_tests OFF CACHE BOOL "")
set(${ExternalProjectName}_build_benchmarks ON CACHE BOOL "")
set(${ExternalProjectName}_install OFF CACHE BOOL "")

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)
