cmake_minimum_required(VERSION 3.9)
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(mba)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
add_subdirectory(python)
add_subdirectory(tests)
add_executable(example example.cpp)
add_test(NAME go_test COMMAND go test WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mba)
