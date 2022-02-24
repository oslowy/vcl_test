cmake_minimum_required(VERSION 3.21)
project(vcl_test)

set(CMAKE_CXX_STANDARD 17)

include_directories(includes/dgemm)
include_directories(includes/vcl)

add_compile_options(
        -O2
)

add_executable(vcl_test main.cpp dgemm_vector.cpp DgemmNonBlocked1AccIJK.cpp includes/dgemm/DgemmNonBlocked1AccIJK.h DgemmNonBlocked1AccKIJ.cpp includes/dgemm/DgemmNonBlocked1AccKIJ.h DgemmBlocked1AccIJK.cpp includes/dgemm/DgemmBlocked1AccIJK.h DgemmBlocked1AccKIJ.cpp includes/dgemm/DgemmBlocked1AccKIJ.h includes/dgemm/dgemm_set.h DgemmVectorNonBlocked.cpp includes/dgemm/DgemmVectorNonBlocked.h DgemmVectorBlocked.cpp includes/dgemm/DgemmVectorBlocked.h DgemmBlocked.cpp)
