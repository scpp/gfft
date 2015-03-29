cmake_minimum_required(VERSION 2.6)

set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(../loki)

#all source files 
add_executable(gfft gfft.cpp)
#add_executable(metapi metapi.cpp)
#add_executable(metasqrt metasqrt.cpp)

# Set default values 
if(NOT DEFINED ${NUM})
set(NUM "16" CACHE STRING "")
endif(NOT DEFINED ${NUM})

if(NOT DEFINED ${FULLOUTPUT})
set(FULLOUTPUT "1" CACHE STRING "")
endif(NOT DEFINED ${FULLOUTPUT})

if(NOT DEFINED ${TYPE})
set(TYPE "DOUBLE" CACHE STRING "")
endif(NOT DEFINED ${TYPE})

if(NOT DEFINED ${PLACE})
set(PLACE "OUT_OF_PLACE" CACHE STRING "")
endif(NOT DEFINED ${PLACE})

if(NOT DEFINED ${NUMTHREADS})
set(NUMTHREADS "1" CACHE STRING "")
endif(NOT DEFINED ${NUMTHREADS})


if(CMAKE_CXX_COMPILER MATCHES "icpc")

### Intel compiler
add_definitions(-openmp -O3 -DGFFTDOC)
add_definitions(-DNUM=${NUM} -DFULLOUTPUT=${FULLOUTPUT} -DTYPE=${TYPE} -DPLACE=${PLACE} -DNUMTHREADS=${NUMTHREADS})

target_link_libraries(gfft iomp5)
#target_link_libraries(metapi iomp5)

elseif(CMAKE_CXX_COMPILER MATCHES "cl.exe")

### MSVC compiler
#add_definitions(/O2 /Ob2 /Oi /Ot /Oy /EHsc /FD /MD /openmp /DGFFTDOC)
add_definitions(/D _USE_MATH_DEFINES)
add_definitions(/D NUM=${NUM} /D FULLOUTPUT=${FULLOUTPUT} /D TYPE=${TYPE} /D PLACE=${PLACE} /D NUMTHREADS=${NUMTHREADS})

elseif(CMAKE_CXX_COMPILER MATCHES "clang")

### Clang compiler
add_definitions(-Wall -O3 -ftemplate-depth=100000 -DFFTW -DQD)
add_definitions(-DNUM=${NUM} -DFULLOUTPUT=${FULLOUTPUT} -DTYPE=${TYPE} -DPLACE=${PLACE} -DNUMTHREADS=${NUMTHREADS})

target_link_libraries(gfft c m stdc++ gomp)
#target_link_libraries(metapi c m stdc++ gomp)

else(CMAKE_CXX_COMPILER MATCHES "icpc")

### GCC compiler
add_definitions(-time -Wall -fopenmp -O0 -g -ftemplate-backtrace-limit=0 -DGFFTDOC)
add_definitions(-DNUM=${NUM} -DFULLOUTPUT=${FULLOUTPUT} -DTYPE=${TYPE} -DPLACE=${PLACE} -DNUMTHREADS=${NUMTHREADS})

target_link_libraries(gfft stdc++ gomp)
# target_link_libraries(metapi c m stdc++ gomp)
# target_link_libraries(metasqrt c m stdc++ gomp)

endif(CMAKE_CXX_COMPILER MATCHES "icpc")

