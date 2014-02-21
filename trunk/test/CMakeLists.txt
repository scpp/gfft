cmake_minimum_required(VERSION 2.6)

project(gfft)


include_directories(../loki)
include_directories(../src)

#add definitions, compiler switches, etc.
#add_definitions(-time -Wall -O0 -g3 -fopenmp -DFFTW -DARPREC)
add_definitions(-time -Wall -O3 -fopenmp -ftemplate-depth=100000 -DFFTW -DQD)
#add_definitions(/O2 /Ob2 /Oi /Ot /Oy /EHsc /FD /MD /openmp /DGFFTDOC)

#list all source files here
add_executable(gfft_performance gfft_performance.cpp)
add_executable(gfft_accuracy gfft_accuracy.cpp)

#need to link to some other libraries ? just add them here
target_link_libraries(gfft_performance gomp)
target_link_libraries(gfft_accuracy gomp fftw3 fftw3l qd)

