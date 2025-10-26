set(ExternalProjectName sr-index)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/sr-index.git
        GIT_TAG feature/cmake
        FIND_PACKAGE_ARGS
)

set(SR-INDEX_ENABLE_TESTS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_TOOLS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_BENCHMARKS ON CACHE BOOL "")
#set(SR-INDEX_INSTALL OFF CACHE BOOL "")

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)
