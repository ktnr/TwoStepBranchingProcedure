# https://cliutils.gitlab.io/modern-cmake/
# for CI: https://raymii.org/s/articles/Github_Actions_cpp_boost_cmake_speedup.html

cmake_minimum_required(VERSION 3.16)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

project(
    tsbp-cli
    VERSION 0.9
    LANGUAGES CXX)
message(STATUS "project: ${PROJECT_NAME}")
message(STATUS "version: ${PROJECT_VERSION}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Projects
add_executable(tsbp 
    "src/tsbp-cli.cpp")

target_include_directories(tsbp PUBLIC
    ${CMAKE_SOURCE_DIR}/external
    ${CMAKE_SOURCE_DIR}/tsbp/include)

if(UNIX)
    target_link_libraries(tsbp PRIVATE tsbp-lib tbb pthread)
else()
    target_link_libraries(tsbp PRIVATE tsbp-lib)
endif(UNIX)