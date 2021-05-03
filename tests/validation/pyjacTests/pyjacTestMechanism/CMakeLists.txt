cmake_minimum_required(VERSION 2.6)

project(pyJac)

set(CMAKE_BUILD_TYPE Release)

enable_language(C)

#Note, this is not the most optimized compiler setup.
set(CMAKE_C_FLAGS "-std=c99 -Ofast -fPIC")

include_directories(out)
include_directories(out/jacobs)

file(GLOB_RECURSE SOURCES "out/*.c")

add_library(c_pyjac_test SHARED ${SOURCES})

install(TARGETS c_pyjac_test DESTINATION .)

