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
download_project(PROJ grammar
        GIT_REPOSITORY https://github.com/duscob/grammar.git
        GIT_TAG master
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})


set(grammar_build_tools ON CACHE BOOL "grammar_build_tools")
set(grammar_build_tests OFF CACHE BOOL "")
set(grammar_build_benchmarks OFF CACHE BOOL "")
set(grammar_install OFF CACHE BOOL "grammar_install")

add_subdirectory(${grammar_SOURCE_DIR} ${grammar_BINARY_DIR})

include_directories("${grammar_SOURCE_DIR}/include")
