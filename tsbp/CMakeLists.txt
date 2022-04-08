# https://cliutils.gitlab.io/modern-cmake/
# for CI: https://raymii.org/s/articles/Github_Actions_cpp_boost_cmake_speedup.html

cmake_minimum_required(VERSION 3.18)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

project(
    tsbp
    VERSION 0.9
    LANGUAGES CXX)
message(STATUS "project: ${PROJECT_NAME}")
message(STATUS "version: ${PROJECT_VERSION}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Boost
#set(Boost_DEBUG 1)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.76 REQUIRED) # REQUIRED COMPONENTS graph dynamic_bitset property_map

message(STATUS "Boost_FOUND: ${Boost_FOUND}")
message(STATUS "Boost version: ${Boost_VERSION}")
message(STATUS "Boost include: ${Boost_INCLUDE_DIR}")

# Projects
add_library(tsbp-lib
    "src/BranchAndBound.cpp"
    "src/Model.cpp"
    "src/OrthogonalPackingSolver2D.cpp"
    "src/Serialization.cpp"
    "src/Preprocess.cpp"
    "include/tsbp/BranchAndBound.h"
    "include/tsbp/Model.h"
    "include/tsbp/OrthogonalPackingSolver2D.h"
    "include/tsbp/Serialization.h"
    "include/tsbp/Preprocess.h"
    )

target_include_directories(tsbp-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/tsbp
    ${CMAKE_SOURCE_DIR}/external
    ${Boost_INCLUDE_DIR})