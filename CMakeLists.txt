cmake_minimum_required(VERSION 3.21)
project(vcl_test)

set(CMAKE_CXX_STANDARD 17)

include_directories(includes/dgemm)
include_directories(includes/vcl)

add_compile_options( -O2 -march=native)

# Increase stack size for large memory allocations
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,40000000")

add_executable(vcl_test main.cpp DgemmScalarNonBlocked1AccIJK.cpp includes/dgemm/DgemmScalarNonBlocked1AccIJK.h DgemmScalarNonBlocked1AccJKI.cpp includes/dgemm/DgemmScalarNonBlocked1AccJKI.h DgemmScalarBlocked1AccIJK.cpp includes/dgemm/DgemmScalarBlocked1AccIJK.h DgemmScalarBlocked1AccJKI.cpp includes/dgemm/DgemmScalarBlocked1AccJKI.h includes/dgemm/DgemmVectorBlocked.h DgemmVectorNonBlocked1AccIJK.cpp includes/dgemm/DgemmVectorNonBlocked1AccIJK.h DgemmVectorNonBlocked1AccJKI.cpp includes/dgemm/DgemmVectorNonBlocked1AccJKI.h DgemmVectorBlocked1AccIJK.cpp includes/dgemm/DgemmVectorBlocked1AccIJK.h DgemmVectorBlocked1AccJKI.cpp includes/dgemm/DgemmVectorBlocked1AccJKI.h DgemmScalarNonBlocked4AccIJK.cpp includes/dgemm/DgemmScalarNonBlocked4AccIJK.h DgemmScalarNonBlocked4AccJKI.cpp includes/dgemm/DgemmScalarNonBlocked4AccJKI.h DgemmScalarBlocked4AccIJK.cpp includes/dgemm/DgemmScalarBlocked4AccIJK.h DgemmScalarBlocked4AccJKI.cpp includes/dgemm/DgemmScalarBlocked4AccJKI.h DgemmReference.cpp includes/dgemm/DgemmReference.h DgemmVectorNonBlocked2AccIJK.cpp includes/dgemm/DgemmVectorNonBlocked2AccIJK.h)
