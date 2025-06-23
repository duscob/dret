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
download_project(PROJ r-index
        GIT_REPOSITORY https://github.com/duscob/r-index.git
        GIT_TAG library
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})


set(rindex_build_tools ON CACHE BOOL "rindex_build_tools")
set(rindex_build_tests OFF CACHE BOOL "")
set(rindex_build_benchmarks OFF CACHE BOOL "")
set(rindex_install OFF CACHE BOOL "rindex_install")

#add_subdirectory(${r-index_SOURCE_DIR} ${r-index_BINARY_DIR})

include_directories("${r-index_SOURCE_DIR}/internal")
