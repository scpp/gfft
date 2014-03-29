cmake_minimum_required(VERSION 2.6)

### If you turn on one of the following compiler definitions
### you have to run cmake twice to take effect

### Compiler
set(CMAKE_CXX_COMPILER icpc)

project(gfft CXX)

#if you don't want the full compiler output, remove the following line
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(db)
