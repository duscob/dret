set(ExternalProjectName sr-index)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/sr-index.git
        GIT_TAG main
        FIND_PACKAGE_ARGS
)

set(SR-INDEX_ENABLE_TESTS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_TOOLS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_BENCHMARKS ON CACHE BOOL "")
#set(SR-INDEX_INSTALL OFF CACHE BOOL "")
set(SR-INDEX_DEFAULT_ALPHABET_WIDTH ${DRET_DEFAULT_ALPHABET_WIDTH} CACHE INTERNAL "")

FetchContent_MakeAvailable(${ExternalProjectName})
FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)

FetchContent_MakeAvailable(json)
include_directories(${json_SOURCE_DIR}/single_include)
