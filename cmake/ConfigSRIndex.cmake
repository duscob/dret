# Adapted from https://github.com/Crascit/DownloadProject/blob/master/CMakeLists.txt
#
# CAVEAT: use DownloadProject.cmake
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (CMAKE_VERSION VERSION_LESS 3.2)
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
else ()
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
endif ()

include(DownloadProject)
download_project(PROJ sr-index
        GIT_REPOSITORY https://github.com/duscob/sr-index.git
        GIT_TAG main
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})


set(SR-INDEX_ENABLE_TESTS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_TOOLS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_BENCHMARKS ON CACHE BOOL "")
#set(SR-INDEX_INSTALL OFF CACHE BOOL "")

add_subdirectory(${sr-index_SOURCE_DIR} ${sr-index_BINARY_DIR})

include_directories("${sr-index_SOURCE_DIR}/include")
