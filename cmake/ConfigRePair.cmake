#include(GNUInstallDirs)
include(ExternalProject)

set(ExternalProjectName repair)

#set(config_flags)  # parameters desired for ./configure of Autotools

#set(${ExternalProjectName}_LIBRARY "${CMAKE_STATIC_LIBRARY_PREFIX}${ExternalProjectName}${CMAKE_STATIC_LIBRARY_SUFFIX}")

set(_src ${PROJECT_BINARY_DIR}/${ExternalProjectName}-prefix/src/${ExternalProjectName}/bal/)

find_program(MAKE_EXECUTABLE NAMES gmake nmake make mingw32-make REQUIRED)

ExternalProject_Add(
        ${ExternalProjectName}
        URL https://users.dcc.uchile.cl/~gnavarro/software/repair.tgz
        DOWNLOAD_EXTRACT_TIMESTAMP true
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE_EXECUTABLE} -j -C ${_src} repair despair irepair idespair
        INSTALL_COMMAND ""
        #        TEST_COMMAND ""
        #        BUILD_BYPRODUCTS <BINARY_DIR>/src/.libs/${${ExternalProjectName}_LIBRARY}
)

ExternalProject_Get_property(${ExternalProjectName} SOURCE_DIR)
set(${ExternalProjectName}_SOURCE_DIR ${SOURCE_DIR})
#include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)

ExternalProject_Get_property(${ExternalProjectName} BINARY_DIR)
set(${ExternalProjectName}_BINARY_DIR ${BINARY_DIR})
