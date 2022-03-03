cmake_minimum_required(VERSION 3.21)
project(vcl_test)

set(CMAKE_CXX_STANDARD 17)

include_directories(includes/dgemm)
include_directories(includes/vcl)

add_compile_options( -O2 -march=native)

# Increase stack size for large memory allocations
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,40000000")

add_executable(vcl_test main.cpp DgemmVector.cpp DgemmNonBlocked1AccIJK.cpp includes/dgemm/DgemmNonBlocked1AccIJK.h DgemmNonBlocked1AccJKI.cpp includes/dgemm/DgemmNonBlocked1AccJKI.h DgemmBlocked1AccIJK.cpp includes/dgemm/DgemmBlocked1AccIJK.h DgemmBlocked1AccJKI.cpp includes/dgemm/DgemmBlocked1AccJKI.h DgemmVectorBlocked.cpp includes/dgemm/DgemmVectorBlocked.h DgemmBlocked.cpp DgemmVectorNonBlockedIJK.cpp includes/dgemm/DgemmVectorNonBlockedIJK.h DgemmVectorNonBlockedJKI.cpp includes/dgemm/DgemmVectorNonBlockedJKI.h DgemmVectorBlockedIJK.cpp includes/dgemm/DgemmVectorBlockedIJK.h DgemmVectorBlockedJKI.cpp includes/dgemm/DgemmVectorBlockedJKI.h DgemmNonBlocked4AccIJK.cpp includes/dgemm/DgemmNonBlocked4AccIJK.h DgemmNonBlocked4AccJKI.cpp includes/dgemm/DgemmNonBlocked4AccJKI.h DgemmBlocked4AccIJK.cpp includes/dgemm/DgemmBlocked4AccIJK.h DgemmBlocked4AccJKI.cpp includes/dgemm/DgemmBlocked4AccJKI.h DgemmReference.cpp includes/dgemm/DgemmReference.h)
