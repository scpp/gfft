
PROJECT(gfft)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -O2 -fopenmp)

#list all source files here
ADD_EXECUTABLE(gfft gfft.cpp)
ADD_EXECUTABLE(bench bench.cpp)

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(gfft gomp)
TARGET_LINK_LIBRARIES(bench gomp)

