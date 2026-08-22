#include(GNUInstallDirs)
include(ExternalProject)

set(ExternalProjectName repair)

#set(config_flags)  # parameters desired for ./configure of Autotools

#set(${ExternalProjectName}_LIBRARY "${CMAKE_STATIC_LIBRARY_PREFIX}${ExternalProjectName}${CMAKE_STATIC_LIBRARY_SUFFIX}")

# Both variants are built, and the choice is made at run time per document array
# (see RunRePair in construct_base.h).
#
# `bal/irepair` indexes the sequence with `int`, so it is correct only while the
# element count fits in one; the largest collections exceed 2^31 entries. The
# distribution ships a `large/` tree (`typedef long long relong`) whose
# `large/bal/` mirrors `bal/`. Given enough <MB> the two emit byte-identical
# .R/.C files, so grammars are interchangeable and no cache is invalidated by
# switching -- which is what makes a per-collection choice safe.
set(_src_bal ${PROJECT_BINARY_DIR}/${ExternalProjectName}-prefix/src/${ExternalProjectName}/bal/)
set(_src_large ${PROJECT_BINARY_DIR}/${ExternalProjectName}-prefix/src/${ExternalProjectName}/large/bal/)

find_program(MAKE_EXECUTABLE NAMES gmake nmake make mingw32-make REQUIRED)

ExternalProject_Add(
        ${ExternalProjectName}
        URL https://users.dcc.uchile.cl/~gnavarro/software/repair.tgz
        DOWNLOAD_EXTRACT_TIMESTAMP true
        # large/bal ships with `gcc -g -m64` and no optimisation at all, which
        # alone makes it 2.2x slower than the -O9 `bal`. With -O9 restored the
        # remaining gap is ~1.15x at 6M ints, ~1.55x at 24M.
        PATCH_COMMAND sed -i "s/gcc -g -m64/gcc -O9 -m64/g" <SOURCE_DIR>/large/bal/makefile
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE_EXECUTABLE} -j -C ${_src_bal} repair despair irepair idespair
              COMMAND ${MAKE_EXECUTABLE} -j -C ${_src_large} repair despair irepair idespair
        INSTALL_COMMAND ""
        #        TEST_COMMAND ""
        #        BUILD_BYPRODUCTS <BINARY_DIR>/src/.libs/${${ExternalProjectName}_LIBRARY}
)

ExternalProject_Get_property(${ExternalProjectName} SOURCE_DIR)
set(${ExternalProjectName}_SOURCE_DIR ${SOURCE_DIR})
#include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)

ExternalProject_Get_property(${ExternalProjectName} BINARY_DIR)
set(${ExternalProjectName}_BINARY_DIR ${BINARY_DIR})
