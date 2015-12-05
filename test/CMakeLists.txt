cmake_minimum_required(VERSION 2.8.7)
project(par)
find_package(PkgConfig REQUIRED)
pkg_search_module(CURL REQUIRED libcurl)

set(CMAKE_C_FLAGS "-std=c99 -O3 -Wall")
set(CMAKE_CXX_FLAGS "-O3 -Wall")

include_directories(
    .
    ..
    ${CURL_INCLUDE_DIRS})

add_executable(
    test_msquares
    test_msquares.c
    lz4.cpp
    lodepng.c
    sds.c
    whereami.c)

target_link_libraries(
    test_msquares
    ${CURL_LIBRARIES})