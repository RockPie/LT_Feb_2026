cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)

project(tb2025_feb_desy)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# Verify conda is not active
if(DEFINED ENV{CONDA_PREFIX})
    message(WARNING "Conda environment detected! This may cause library conflicts.")
    message(WARNING "Consider deactivating conda before running cmake.")
endif()

# Set ROOT installation path explicitly
set(ROOT_DIR "/home/shihai/sw/root/root_install") # Adjust to your ROOT installation path
list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})

# Configure RPATH to only use custom ROOT installation
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_INSTALL_RPATH "${ROOT_DIR}/lib")

# Add linker flags to prefer custom ROOT
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,${ROOT_DIR}/lib")

# Add optimization flags and OpenMP
find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP NOT found! Parallelization will NOT be enabled.")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

# Find ROOT with correct version
find_package(ROOT REQUIRED COMPONENTS
    Core 
    RIO
    Tree
    Gpad    
    Hist
    Thread
    MathCore
    ROOTDataFrame 
)

include_directories(${ROOT_INCLUDE_DIRS})

# Ensure FetchContent is available
include(FetchContent)

# Fetch argparse
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG        master  # Consider specifying a stable release tag
)
FetchContent_MakeAvailable(argparse)

# Fetch easylogging++
FetchContent_Declare(
    easyloggingpp
    GIT_REPOSITORY https://github.com/muflihun/easyloggingpp.git
    GIT_TAG        v9.97.0
)
FetchContent_MakeAvailable(easyloggingpp)

# Fetch nlohmann/json
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3  # Latest stable release
)
FetchContent_MakeAvailable(json)

# Fetch CSV for Modern C++
FetchContent_Declare(
    csv_parser
    GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser.git
    GIT_TAG        2.1.3  # latest stable release, adjust as needed
)
FetchContent_MakeAvailable(csv_parser)

# Include easylogging++ headers
include_directories(${easyloggingpp_SOURCE_DIR}/src)

# Compile easylogging++.cc as a static library
add_library(easyloggingpp STATIC ${easyloggingpp_SOURCE_DIR}/src/easylogging++.cc)

# Add common library
add_library(H2GCROC_Common STATIC
    src/H2GCROC_Common.cxx
)

target_include_directories(H2GCROC_Common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${easyloggingpp_SOURCE_DIR}/src
    ${ROOT_INCLUDE_DIRS}
    ${json_SOURCE_DIR}/include
    ${csv_parser_SOURCE_DIR}/single_include
)

target_link_libraries(H2GCROC_Common PUBLIC
    easyloggingpp
    argparse
    nlohmann_json::nlohmann_json  # Link JSON library
    ${ROOT_LIBRARIES}
    ROOT::ROOTDataFrame
    OpenMP::OpenMP_CXX
)

# Add lib
add_library(H2GCROC_Lib STATIC
    src/H2GCROC_Lib.cxx
)

target_include_directories(H2GCROC_Lib PUBLIC
     ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${easyloggingpp_SOURCE_DIR}/src
    ${ROOT_INCLUDE_DIRS}
    ${json_SOURCE_DIR}/include
    ${csv_parser_SOURCE_DIR}/single_include
)

target_link_libraries(H2GCROC_Lib PUBLIC
    H2GCROC_Common
    easyloggingpp
    argparse
    nlohmann_json::nlohmann_json  # Link JSON library
    ${ROOT_LIBRARIES}
    ROOT::ROOTDataFrame
    OpenMP::OpenMP_CXX
)


# Find all script executables
file(GLOB scripts ${CMAKE_CURRENT_SOURCE_DIR}/scripts/*.cxx)
foreach(script ${scripts})
    get_filename_component(script_name ${script} NAME_WE)
    add_executable(${script_name} ${script})

    target_include_directories(${script_name} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${easyloggingpp_SOURCE_DIR}/src
        ${json_SOURCE_DIR}/include
        ${csv_parser_SOURCE_DIR}/single_include
    )

    target_link_libraries(${script_name}
        H2GCROC_Common
        H2GCROC_Lib
        easyloggingpp
        argparse
        nlohmann_json::nlohmann_json
        ${ROOT_LIBRARIES}
        ROOT::ROOTDataFrame
        OpenMP::OpenMP_CXX
    )
endforeach(script)