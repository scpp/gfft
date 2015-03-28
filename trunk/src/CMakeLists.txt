cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(../loki)

#all source files 
#add_executable(gfft gfft.cpp)
#add_executable(gfft_inp gfft_inp.cpp)
#add_executable(gfft_oop gfft_oop.cpp)
add_executable(cgfft_inp cgfft_inp.cpp)
add_executable(cgfft_oop cgfft_oop.cpp)
#add_executable(metapi metapi.cpp)
#add_executable(metasqrt metasqrt.cpp)

#add_executable(metaf metaf.cpp)

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

elseif(CMAKE_CXX_COMPILER MATCHES "clang")

### Clang compiler
add_definitions(-Wall -O3 -ftemplate-depth=100000 -DFFTW -DQD)

target_link_libraries(gfft_inp c m stdc++ gomp)
target_link_libraries(gfft_oop c m stdc++ gomp)
target_link_libraries(cgfft_inp c m stdc++ gomp)
target_link_libraries(cgfft_oop c m stdc++ gomp)
target_link_libraries(metapi c m stdc++ gomp)

else(CMAKE_CXX_COMPILER MATCHES "icpc")

### GCC compiler
add_definitions(-time -Wall -DNUM=${MYNUM} -DAUTO=${MYAUTO} -fopenmp -O0 -g -ftemplate-backtrace-limit=0 -DGFFTDOC)

#target_link_libraries(gfft stdc++ gomp)
#target_link_libraries(gfft_inp stdc++ gomp)
#target_link_libraries(gfft_oop stdc++ gomp)
target_link_libraries(cgfft_inp stdc++ gomp)
target_link_libraries(cgfft_oop stdc++ gomp)
# target_link_libraries(metapi c m stdc++ gomp)
# target_link_libraries(metasqrt c m stdc++ gomp)

#target_link_libraries(metaf c m stdc++ gomp)

endif(CMAKE_CXX_COMPILER MATCHES "icpc")

