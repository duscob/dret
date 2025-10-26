set(ExternalProjectName json)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.12.0
        #        DOWNLOAD_EXTRACT_TIMESTAMP true
        #        URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
)

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/single_include)
