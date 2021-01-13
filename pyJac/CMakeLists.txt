cmake_minimum_required(VERSION 2.6)

project(pyJac)

set(CMAKE_BUILD_TYPE Release)

enable_language(C)

#GNU on DESKTOP
set(CMAKE_C_FLAGS "-std=c99 -O3 -mtune=native -fPIC")

#GNU on SISU
#set(CMAKE_C_FLAGS "-std=c99 -fPIC")
#set(CMAKE_C_FLAGS "-std=c99 -O3 -ffast-math -funroll-loops -march=haswell -mtune=haswell -fPIC")
#set(CMAKE_C_FLAGS "-std=c99 -O3 -ffast-math -funroll-loops -march=haswell -mtune=haswell -fPIC -mavx2 -mfma")
#set(CMAKE_C_FLAGS "-std=c99 -O3 -ffast-math -funroll-loops -march=haswell -mtune=haswell -fPIC -mavx2 -mfma")


include_directories(out)
include_directories(out/jacobs)


file(GLOB_RECURSE SOURCES "out/*.c")

add_library(c_pyjac SHARED ${SOURCES})

install(TARGETS c_pyjac DESTINATION .)

