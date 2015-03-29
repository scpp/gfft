cmake_minimum_required(VERSION 2.6)

include_directories(../loki)
include_directories(../src)

#list all source files here
add_executable(gfft_performance gfft_performance.cpp)
#add_executable(gfft_accuracy gfft_accuracy.cpp)

### Compiler
set(CMAKE_CXX_COMPILER g++)


# Set default values 
if(NOT DEFINED ${PNUM})
set(PNUM "2" CACHE STRING "")
endif(NOT DEFINED ${PNUM})

if(NOT DEFINED ${FULLOUTPUT})
set(FULLOUTPUT "1" CACHE STRING "")
endif(NOT DEFINED ${FULLOUTPUT})

if(NOT DEFINED ${PMIN})
set(PMIN "2" CACHE STRING "")
endif(NOT DEFINED ${PMIN})

if(NOT DEFINED ${PMAX})
set(PMAX "4" CACHE STRING "")
endif(NOT DEFINED ${PMAX})

if(NOT DEFINED ${TYPE})
set(TYPE "DOUBLE" CACHE STRING "")
endif(NOT DEFINED ${TYPE})

if(NOT DEFINED ${PLACE})
set(PLACE "OUT_OF_PLACE" CACHE STRING "")
endif(NOT DEFINED ${PLACE})


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
add_definitions(-DPNUM=${PNUM} -DFULLOUTPUT=${FULLOUTPUT} -DPMIN=${PMIN} -DPMAX=${PMAX} -DTYPE=${TYPE} -DPLACE=${PLACE})

target_link_libraries(gfft_performance stdc++ qd gomp)
#target_link_libraries(gfft_accuracy stdc++ gomp qd fftw3 fftw3l)

endif(CMAKE_CXX_COMPILER MATCHES "icpc")
