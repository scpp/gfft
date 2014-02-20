cmake_minimum_required(VERSION 2.6)

project(gfft)

#if you don't want the full compiler output, remove the following line
set(CMAKE_VERBOSE_MAKEFILE ON)

### If you turn on one of the following compiler definitions
### you have to run cmake twice to take effect
### GNU GCC compiler
#set(CMAKE_CXX_COMPILER g++)
#set(CMAKE_CC_COMPILER gcc)

### Intel compiler
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_CC_COMPILER icc)

#add_subdirectory(src)
add_subdirectory(test)
