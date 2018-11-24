find_package(OpenMP)

add_executable(test_mba test_mba.cpp)
target_link_libraries(test_mba OpenMP::OpenMP_CXX)
add_test(test_mba test_mba)

