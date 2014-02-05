cmake_minimum_required(VERSION 2.6)

project(gfft)


include_directories(../loki)

### GCC compiler
add_definitions(-time -Wall -O3 -fopenmp -ftemplate-depth=10000 -DGFFTDOC)
#add_definitions(-time -Wall -O0 -g3 -fopenmp -ftemplate-depth=10000 -DGFFTDOC)

### Intel compiler
#add_definitions(-O3 -openmp -DGFFTDOC)

### MSVC compiler
#add_definitions(/O2 /Ob2 /Oi /Ot /Oy /EHsc /FD /MD /openmp /DGFFTDOC)

#list all source files here
#add_executable(gfft gfft.cpp)
add_executable(gfft_inp gfft_inp.cpp)
#add_executable(gfft_oop gfft_oop.cpp)
#add_executable(cgfft_inp cgfft_inp.cpp)
#add_executable(cgfft_oop cgfft_oop.cpp)
#add_executable(metapi metapi.cpp)

### Use following lines for GCC compiler
#target_link_libraries(gfft gomp)
target_link_libraries(gfft_inp gomp)
#target_link_libraries(gfft_oop gomp)
#target_link_libraries(cgfft_inp gomp)
#target_link_libraries(cgfft_oop gomp)

### Use following lines for Intel compiler
#target_link_libraries(gfft iomp5 pthread)
#target_link_libraries(cgfft iomp5 pthread)

