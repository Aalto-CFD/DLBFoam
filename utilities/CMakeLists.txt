cmake_minimum_required(VERSION 2.6)

project(pyJac)

set(CMAKE_BUILD_TYPE Release)

enable_language(C)

set(CMAKE_C_FLAGS "-std=c99 -Ofast -march=native -fPIC")

include_directories(out)
include_directories(out/jacobs)


file(GLOB_RECURSE SOURCES "out/*.c")

add_library(c_pyjac SHARED ${SOURCES})

install(TARGETS c_pyjac DESTINATION .)

