set(ExternalProjectName gflags)

include(FetchContent)
# Deliberately no FIND_PACKAGE_ARGS: always build from source, which exports
# gflags::gflags. An installed gflags 2.2.x instead exports a bare imported
# gflags target, so allowing find_package here makes the target name vary by host.
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/gflags/gflags.git
        GIT_TAG v2.3.0
)

set(GFLAGS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(GFLAGS_BUILD_PACKAGING OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(${ExternalProjectName})
