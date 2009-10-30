cmake_minimum_required(VERSION 2.6)

project(gfft)

#if you don't want the full compiler output, remove the following line
set(CMAKE_VERBOSE_MAKEFILE ON)

#include_directories()

#add definitions, compiler switches, etc.
add_definitions(-time -Wall -O2 -fopenmp -DGFFTDOC -DFFTW)

#list all source files here
add_executable(gfft gfft.cpp)
#add_executable(cgfft cgfft.cpp)
#add_executable(bench bench.cpp)
#add_executable(ffttest ffttest.cpp)
#add_executable(check check.cpp)

#need to link to some other libraries ? just add them here
target_link_libraries(gfft gomp)
#target_link_libraries(cgfft gomp)
#target_link_libraries(bench gomp)

