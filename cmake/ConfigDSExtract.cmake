include(ExternalProject)

set(ExternalProjectName dsextract)

# Cecilia Hernandez's web-graph dense-subgraph (biclique) extractor used by
# drl's PDL-BC. Source vendored at external/dsextract/ — see NOTES.md
# there for provenance and the boost-removal patch.
#
# We don't fetch from a URL because no upstream URL is known; the
# tarball lives in drl's tree (drl/external/doclist-env/software/
# dsextract.tgz) and we ship the already-patched source directly.

set(_vendored ${PROJECT_SOURCE_DIR}/external/dsextract)
set(_src ${PROJECT_BINARY_DIR}/${ExternalProjectName}-prefix/src/${ExternalProjectName}/)

find_program(MAKE_EXECUTABLE NAMES gmake nmake make mingw32-make REQUIRED)

ExternalProject_Add(
        ${ExternalProjectName}
        SOURCE_DIR ${_vendored}
        DOWNLOAD_COMMAND ""
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE FALSE
        # The upstream Makefile builds in-tree; copy the vendored sources
        # into the per-build src/ dir on every configure so the build
        # output stays out of the source tree.
        BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${_vendored} ${_src}
                COMMAND ${MAKE_EXECUTABLE} -j -C ${_src} vnmextract
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS ${_src}/vnmextract
)

ExternalProject_Get_property(${ExternalProjectName} SOURCE_DIR)
set(${ExternalProjectName}_SOURCE_DIR ${SOURCE_DIR})

# vnmextract ends up in the per-build copy of the source dir, not the
# vendored dir, because the Makefile builds in place.
set(${ExternalProjectName}_VNMEXTRACT ${_src}/vnmextract)
