if(CMAKE_VERSION VERSION_LESS "3.5")
	cmake_minimum_required(VERSION 3.1)
else()
	# Ubuntu 26.04 ships CMake > 4, which requires project compatibility >= 3.5.
	cmake_minimum_required(VERSION 3.5)
endif()

project(pyJac)

set(CMAKE_BUILD_TYPE Release)

enable_language(C)

# Note, this is not the most optimized compiler setup.
set(CMAKE_C_FLAGS "-std=c99 -Ofast -fPIC")

include_directories(src)
include_directories(src/jacobs)

file(GLOB_RECURSE SOURCES "src/*.c")

add_library(c_pyjac_test SHARED ${SOURCES})

install(TARGETS c_pyjac_test DESTINATION .)
