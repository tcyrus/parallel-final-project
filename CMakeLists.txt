cmake_minimum_required(VERSION 3.10)
project(SEIR C)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_COMPILER mpicc)
#set(CMAKE_C_COMPILER mpixlc)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread -Wall -Wextra")

add_executable(SEIR SEIR.c cell_state.h person.h board.h)
#add_executable(SEIR.xl SEIR.c)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

target_link_libraries(SEIR Threads::Threads)
target_link_libraries(SEIR ${MPI_C_LIBRARIES})
