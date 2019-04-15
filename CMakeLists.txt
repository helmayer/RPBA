# RPBA - Robust Parallel Bundle Adjustment
#
# File CMakeLists.txt
#
#
#
# Copyright 2019 Helmut Mayer, Bundeswehr University Munich, Germany, Helmut.Mayer@unibw.de
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


cmake_minimum_required(VERSION 3.5)
if(WIN32)
	project(rpba CXX)
else()
	project(rpba)
endif(WIN32)


# RPBA files
file(GLOB rpba_sources rpba.cpp rpbacore.cpp rbacore.cpp)


set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "limited configs" FORCE)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build mode (Values: Release, Debug)")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)


# ------------------------------------
# METIS
# ------------------------------------
if(CMAKE_COMPILER_IS_GNUCXX)
	set(METIS_LIBS libmetis.so)
else()
	set(METIS_LIBS metis.lib)
endif()
find_path(METIS_INCLUDE_DIR metis.h $ENV{METIS_ROOT}/include /usr/include)
find_path(METIS_LIB_DIR ${METIS_LIBS} $ENV{METIS_ROOT}/lib /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu)
include_directories(${METIS_INCLUDE_DIR})
link_directories(${METIS_LIB_DIR})
set(HPE_TARGET_LIBS ${HPE_TARGET_LIBS} ${METIS_LIBS})


if(UNIX)
	option(BUILD_WITH_TCMALLOC "Use TCMalloc" OFF)
	option(BUILD_WITH_PROFILER "Use TCMalloc with profiler" OFF)

	# TCMalloc
	if(BUILD_WITH_TCMALLOC)
		if(BUILD_WITH_PROFILER)
			set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} tcmalloc_and_profiler)
		else()
			set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} tcmalloc)
		endif()
	endif()

	# TCMalloc with profiler
	if(BUILD_WITH_PROFILER)
		if(BUILD_WITH_TCMALLOC)
		else()
			set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} tcmalloc_and_profiler)
		endif()
	endif()
endif()


# set OpenMP flags
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Check for the 64-bit compiler
if( CMAKE_SIZEOF_VOID_P LESS 8 )
	message(SEND_ERROR "Only 64 bit build mode is supported.")
endif()


include_directories(.)
if(CMAKE_COMPILER_IS_GNUCXX)
#-Wconversion
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unknown-pragmas -Wextra -Wunreachable-code -std=c++14 -m64")
	set(CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-result -Wno-unused-parameter -Wno-unused-local-typedefs -pipe -march=native -pedantic -O3 -DNDEBUG -DEIGEN_NO_DEBUG")
	set(CMAKE_CXX_FLAGS_DEBUG "-Wconversion -pipe -march=native -pedantic -g -DHPE_RUNTIME_CHECKS")
else()
	add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4251 /wd4512")
	set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} "/MP /Ox /DNDEBUG /DEIGEN_NO_DEBUG")
	set(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} "/MP /Ob1")
endif()
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
set(CMAKE_C_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})

if(CMAKE_COMPILER_IS_GNUCXX)
  set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} metis.so)
endif()
      
add_executable(rpba ${rpba_sources})
target_link_libraries(rpba ${rpba_TARGET_LIBS})
