cmake_minimum_required(VERSION 2.6)

project(gfft)

#if you don't want the full compiler output, remove the following line
set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(../src)

#add definitions, compiler switches, etc.
add_definitions(-time -Wall -O3 -fopenmp -DFFTW)

#list all source files here
add_executable(bench bench.cpp)
#add_executable(check check.cpp)

#need to link to some other libraries ? just add them here
target_link_libraries(bench gomp)
#target_link_libraries(check gomp fftw3)

