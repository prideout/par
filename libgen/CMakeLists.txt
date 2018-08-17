cmake_minimum_required(VERSION 2.8.7)
project(libpar)
#find_package(PkgConfig REQUIRED)
#pkg_search_module(CURL REQUIRED libcurl)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-std=c11 -Wall")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wvla -Wall")

include_directories(..)

add_library(shapes MODULE shapes.c)
add_library(msquares MODULE msquares.c)

