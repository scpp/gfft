cmake_minimum_required(VERSION 2.6)

include_directories(../loki)
include_directories(../src)

#list all source files here
add_executable(gfft_performance gfft_performance.cpp)
add_executable(gfft_accuracy gfft_accuracy.cpp)

### Compiler
set(CMAKE_CXX_COMPILER g++)


if(CMAKE_CXX_COMPILER MATCHES "icpc")

### Intel compiler
add_definitions(-O3 -openmp -fno-alias -unroll -DFFTW -DQD)

target_link_libraries(gfft_performance iomp5)
target_link_libraries(gfft_accuracy iomp5 fftw3 fftw3l qd)

elseif(CMAKE_CXX_COMPILER MATCHES "cl.exe")

### MSVC compiler
add_definitions(/O2 /Ob2 /Oi /Ot /Oy /EHsc /FD /MD /openmp /DGFFTDOC)

elseif(CMAKE_CXX_COMPILER MATCHES "clang")

### Clang compiler
add_definitions(-Wall -O3 -ftemplate-depth=100000 -DFFTW -DQD)

target_link_libraries(gfft_performance c m stdc++ qd gomp)
target_link_libraries(gfft_accuracy c m stdc++ gomp qd fftw3 fftw3l)

else(CMAKE_CXX_COMPILER MATCHES "icpc")

### GCC compiler
add_definitions(-time -Wall -O3 -fopenmp -ftemplate-depth=10000 -DFFTW -DQD)

target_link_libraries(gfft_performance stdc++ qd gomp)
target_link_libraries(gfft_accuracy stdc++ gomp qd fftw3 fftw3l)

endif(CMAKE_CXX_COMPILER MATCHES "icpc")
