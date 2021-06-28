cmake_minimum_required(VERSION 2.6)

project(pyJac)

set(CMAKE_BUILD_TYPE Release)

enable_language(C)

#Note, this is not the most optimized compiler setup.
set(CMAKE_C_FLAGS "-std=c99 -Ofast -fPIC")

include_directories(src)
include_directories(src/jacobs)

file(GLOB_RECURSE SOURCES "src/*.c")

add_library(c_pyjac_test SHARED ${SOURCES})

install(TARGETS c_pyjac_test DESTINATION .)
