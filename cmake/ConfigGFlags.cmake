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

set(ExternalProjectName gflags)

include(DownloadProject)
download_project(PROJ ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/gflags/gflags.git
        GIT_TAG v2.2.2
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})

set(GFLAGS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(GFLAGS_BUILD_PACKAGING OFF CACHE BOOL "" FORCE)
set(GFLAGS_BUILD_PACKAGING OFF CACHE BOOL "" FORCE)

add_subdirectory(${${ExternalProjectName}_SOURCE_DIR} ${${ExternalProjectName}_BINARY_DIR} EXCLUDE_FROM_ALL)

include_directories("${${ExternalProjectName}_SOURCE_DIR}/include")
