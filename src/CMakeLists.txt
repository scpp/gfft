cmake_minimum_required(VERSION 2.6)

project(gfft)

#if you don't want the full compiler output, remove the following line
set(CMAKE_VERBOSE_MAKEFILE ON)

#include_directories()

#add definitions, compiler switches, etc.
#add_definitions(-time -Wall -O3 -fopenmp -DGFFTDOC)
add_definitions(-O3 -openmp -DGFFTDOC)

#list all source files here
add_executable(gfft gfft.cpp)
add_executable(cgfft cgfft.cpp)

#need to link to some other libraries ? just add them here
target_link_libraries(gfft gomp)
target_link_libraries(cgfft gomp)
#target_link_libraries(gfft iomp5 pthread)
#target_link_libraries(cgfft iomp5 pthread)

