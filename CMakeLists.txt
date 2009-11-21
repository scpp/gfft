cmake_minimum_required(VERSION 2.6)

project(gfft)

#if you don't want the full compiler output, remove the following line
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CC_COMPILER gcc)


add_subdirectory(src)
#add_subdirectory(test)
