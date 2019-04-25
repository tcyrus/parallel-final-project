cmake_minimum_required(VERSION 3.10)
project(SEIR C)

#find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER gcc)
#set(CMAKE_C_COMPILER mpicc)
#set(CMAKE_C_COMPILER mpixlc)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread -Wall -Wextra -O5")

add_executable(SEIR SEIR.c)
#add_executable(SEIR.xl SEIR.c)

#include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#target_link_libraries(SEIR ${MPI_C_LIBRARIES})
