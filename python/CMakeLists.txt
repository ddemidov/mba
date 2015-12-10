find_package(PythonInterp)
find_package(PythonLibs ${PYTHON_VERSION_STRING})
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(NumPy QUIET)
if (NUMPY_FOUND)
    include_directories(${CMAKE_SOURCE_DIR})
    include_directories(${NUMPY_INCLUDE_DIRS})

    set(BOOST_PYTHON_COMPONENT "python" CACHE STRING "Boost python component (python/python3/...)")

    find_package(Boost COMPONENTS ${BOOST_PYTHON_COMPONENT})

    add_library(pymba SHARED pymba.cpp)
    set_target_properties(pymba PROPERTIES PREFIX "")
    target_link_Libraries(pymba ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

