cmake_minimum_required(VERSION 3.18)

project(
    tsbp-tests
    VERSION 0.9
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/external/doctest/scripts/cmake)

include(doctest)

add_executable(${PROJECT_NAME} 
    "main.cpp"
    "Preprocess.cpp"
    "PlacementLMAO.cpp")

target_include_directories(tsbp-tests PUBLIC
    ${CMAKE_SOURCE_DIR}/external
    ${CMAKE_SOURCE_DIR}/external/doctest/doctest
    ${CMAKE_SOURCE_DIR}/tsbp/include)

target_link_libraries(tsbp-tests PRIVATE tsbp-lib)

# https://github.com/doctest/doctest/blob/master/doc/markdown/build-systems.md
doctest_discover_tests(${PROJECT_NAME})
