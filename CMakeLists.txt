cmake_minimum_required(VERSION 3.17)
project(prefix_scan C)

set(CMAKE_C_STANDARD 11)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(prefix_scan main.c)
target_link_libraries(prefix_scan ${MPI_LIBRARIES})
if(MPI_COMPILE_FLAGS)
    set_target_properties(prefix_scan PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(prefix_scan PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()