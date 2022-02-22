cmake_minimum_required(VERSION 3.21)
project(vcl_test)

set(CMAKE_CXX_STANDARD 17)

include_directories(includes/dgemm)
include_directories(includes/vcl)

add_executable(vcl_test main.cpp includes/dgemm/dgemm_naive.cpp includes/dgemm/dgemm_blocked_vector.cpp includes/dgemm/dgemm_blocked.cpp)
