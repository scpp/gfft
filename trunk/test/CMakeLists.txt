cmake_minimum_required(VERSION 2.6)

include_directories(../loki)
include_directories(../src)

#list all source files here
add_executable(gfft_performance gfft_performance.cpp)
add_executable(gfft_accuracy gfft_accuracy.cpp)

#need to link to some other libraries ? just add them here

if(CMAKE_CXX_COMPILER MATCHES "icpc")

### Intel compiler
add_definitions(-openmp -O3 -DFFTW -DQD)

target_link_libraries(gfft_performance iomp5)
target_link_libraries(gfft_accuracy iomp5 fftw3 fftw3l qd)

elseif(CMAKE_CXX_COMPILER MATCHES "cl.exe")

### MSVC compiler
add_definitions(/O2 /Ob2 /Oi /Ot /Oy /EHsc /FD /MD /openmp /DGFFTDOC)

else(CMAKE_CXX_COMPILER MATCHES "icpc")

### GCC compiler 
add_definitions(-time -Wall -O3 -fopenmp -ftemplate-depth=100000 -DFFTW -DQD)

target_link_libraries(gfft_performance gomp)
target_link_libraries(gfft_accuracy gomp fftw3 fftw3l qd)

endif(CMAKE_CXX_COMPILER MATCHES "icpc")
