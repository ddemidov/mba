set(MBA_PYTHON_VERSION "" CACHE STRING "Python version to use for compiling the MBA library")

# Try to autodetect Python (can be overridden manually if needed)
set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
find_package(PythonLibs ${MBA_PYTHON_VERSION} REQUIRED)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Enable C++11 mode on C++ / Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  # Enable link time optimization and set the default symbol
  # visibility to hidden (very important to obtain small binaries)
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    # Default symbol visibility
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

    # Check for Link Time Optimization support
    CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
    if (HAS_LTO_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    endif()
  endif()
endif()

# Include path for Python header files
include_directories(${PYTHON_INCLUDE_DIR})

# Include path for pybind11 header files -- this may need to be changed depending on your setup
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include)

# Create the binding library
add_library(mba SHARED pymba.cpp)

# Don't add a 'lib' prefix to the shared library
set_target_properties(mba PROPERTIES PREFIX "")

if (WIN32)
  if (MSVC)
    # Enforce size-based optimization and link time code generation
    # on MSVC (~30% smaller binaries in experiments). /bigobj is needed
    # for bigger binding projects due to the limit to 64k addressable sections
    # /MP enables multithreaded builds (relevant when there are many files).
    set_target_properties(mba PROPERTIES COMPILE_FLAGS "/Os /GL /MP /bigobj")
    set_target_properties(mba PROPERTIES LINK_FLAGS "/LTCG")
  endif()

  # .PYD file extension on Windows
  set_target_properties(mba PROPERTIES SUFFIX ".pyd")

  # Link against the Python shared library
  target_link_libraries(mba ${PYTHON_LIBRARY})
elseif (UNIX)
  # It's quite common to have multiple copies of the same Python version
  # installed on one's system. E.g.: one copy from the OS and another copy
  # that's statically linked into an application like Blender or Maya.
  # If we link our plugin library against the OS Python here and import it
  # into Blender or Maya later on, this will cause segfaults when multiple
  # conflicting Python instances are active at the same time.

  # Windows is not affected by this issue since it handles DLL imports
  # differently. The solution for Linux and Mac OS is simple: we just don't
  # link against the Python library. The resulting shared library will have
  # missing symbols, but that's perfectly fine -- they will be resolved at
  # import time.

  # .SO file extension on Linux/Mac OS
  set_target_properties(mba PROPERTIES SUFFIX ".so")

  # Strip unnecessary sections of the binary on Linux/Mac OS
  if(APPLE)
    set_target_properties(mba PROPERTIES MACOSX_RPATH ".")
    set_target_properties(mba PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -dead_strip")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET mba POST_BUILD COMMAND strip -u -r ${PROJECT_BINARY_DIR}/mba.so)
    endif()
  else()
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET mba POST_BUILD COMMAND strip ${PROJECT_BINARY_DIR}/mba.so)
    endif()
  endif()
endif()

