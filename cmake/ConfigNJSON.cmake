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
download_project(PROJ json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})


#add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR})

#include_directories("${json_SOURCE_DIR}/include")
include_directories("${json_SOURCE_DIR}/single_include")

