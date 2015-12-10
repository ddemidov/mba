include_directories(${CMAKE_SOURCE_DIR})

add_definitions(-DMBA_VERBOSE)

add_executable(1d        1d.cpp)
add_executable(2d        2d.cpp)
add_executable(serialize serialize.cpp)
