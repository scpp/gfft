cmake_minimum_required(VERSION 2.6)

set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(../loki)

#all source files 
add_executable(gfft_inp gfft_inp.cpp)
add_executable(gfft_oop gfft_oop.cpp)
add_executable(cgfft_inp cgfft_inp.cpp)
add_executable(cgfft_oop cgfft_oop.cpp)
add_executable(metapi metapi.cpp)

if(CMAKE_CXX_COMPILER MATCHES "icpc")

### Intel compiler
add_definitions(-openmp -O3 -DGFFTDOC)

target_link_libraries(gfft_inp iomp5)
target_link_libraries(gfft_oop iomp5)
target_link_libraries(cgfft_inp iomp5)
target_link_libraries(cgfft_oop iomp5)
target_link_libraries(metapi iomp5)

elseif(CMAKE_CXX_COMPILER MATCHES "cl.exe")

### MSVC compiler
add_definitions(/O2 /Ob2 /Oi /Ot /Oy /EHsc /FD /MD /openmp /DGFFTDOC)

else(CMAKE_CXX_COMPILER MATCHES "icpc")

### GCC compiler 
add_definitions(-time -Wall -fopenmp -O3 -DGFFTDOC)

target_link_libraries(gfft_inp gomp)
target_link_libraries(gfft_oop gomp)
target_link_libraries(cgfft_inp gomp)
target_link_libraries(cgfft_oop gomp)

endif(CMAKE_CXX_COMPILER MATCHES "icpc")

