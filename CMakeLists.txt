cmake_minimum_required(VERSION 2.8)
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
endif()

project(mba)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
add_subdirectory(python)
add_subdirectory(tests)
