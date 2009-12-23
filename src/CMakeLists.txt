cmake_minimum_required(VERSION 2.6)

project(gfft)


#include_directories()

### GCC compiler
add_definitions(-time -Wall -O3 -fopenmp -DGFFTDOC)

### Intel compiler
#add_definitions(-O3 -openmp -DGFFTDOC)

#list all source files here
add_executable(gfft gfft.cpp)
add_executable(cgfft cgfft.cpp)

### Use following lines for GCC compiler
target_link_libraries(gfft gomp)
target_link_libraries(cgfft gomp)

### Use following lines for Intel compiler
#target_link_libraries(gfft iomp5 pthread)
#target_link_libraries(cgfft iomp5 pthread)

