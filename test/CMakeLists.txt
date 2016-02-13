cmake_minimum_required(VERSION 2.8.7)
project(par)
find_package(PkgConfig REQUIRED)
pkg_search_module(CURL REQUIRED libcurl)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-std=c11 -Wall")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wvla -Wall")

include_directories(.. . ${CURL_INCLUDE_DIRS})

add_executable(
    test_msquares
    test_msquares.c
    lz4.cpp
    lodepng.c
    sds.c
    whereami.c)
target_link_libraries(test_msquares ${CURL_LIBRARIES})

add_executable(
    test_bluenoise
    test_bluenoise.c
    lz4.cpp
    lodepng.c
    sds.c
    whereami.c)
target_link_libraries(test_bluenoise ${CURL_LIBRARIES})

add_executable(
    test_cpp
    test_cpp.cpp
    lz4.cpp
    lodepng.c
    sds.c
    whereami.c)
target_link_libraries(test_cpp ${CURL_LIBRARIES})

add_executable(
    test_linkage
    test_linkage.cpp
    test_linkage.c
    lz4.cpp
    lodepng.c
    sds.c
    whereami.c)
target_link_libraries(test_linkage ${CURL_LIBRARIES})

add_executable(
    test_shapes
    test_shapes.c
    console-colors.c)
target_link_libraries(test_shapes m)

add_executable(
    test_bubbles
    test_bubbles.c
    console-colors.c)
target_link_libraries(test_bubbles m)

add_executable(
    test_filecache
    test_filecache.c
    console-colors.c)

add_executable(
    test_filecache_lz4
    test_filecache_lz4.c
    console-colors.c)

add_executable(
    test_sprune
    test_sprune.c
    console-colors.c)
