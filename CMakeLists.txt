cmake_minimum_required(VERSION 3.16)

project(TwoStepBranchingProcedure VERSION 0.9 LANGUAGES CXX)

# Build Dependencies
add_subdirectory(tsbp)
add_subdirectory(tsbp-cli)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(tests)
endif()