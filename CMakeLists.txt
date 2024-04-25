cmake_minimum_required(VERSION 3.27)
project(FastDCD C)

set(CMAKE_C_STANDARD 23)

add_executable(FastDCD main.c
        matrix_operations/Standard_operations.h
        matrix_operations/Standard_operations.c
        matrix_operations/DCBF_operations.h
        matrix_operations/DCBF_operations.c
        echo_canceling.c
        echo_canceling.h
        operations_unit_tests.c
        operations_unit_tests.h
        standard_RLSDCD/standard_RLSDCD.c
        standard_RLSDCD/standard_RLSDCD.h
        DCBF_RLSDCD/DCBF_RLSDCD.c
        DCBF_RLSDCD/DCBF_RLSDCD.h
)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(FastDCD ${GSL_LIBRARIES})