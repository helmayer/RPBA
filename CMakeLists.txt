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
	project(cls CXX)
else()
	project(cls)
endif(WIN32)




# CLS files
file(GLOB cls_sources cls.cpp ../hpen/cpls.cpp ../hpen/clsbase.cpp ../hpen/Camera.cpp ../hpen/cmath.cpp ../hpen/Params.cpp ../hpen/System.cpp ../hpen/Util.cpp)
#foreach(f ${cls_sources})
#	if((${f} MATCHES "- Kopie") OR (${f} MATCHES "sim.cpp"))
#		#message(STATUS "Exclude [" ${f} "] from project.")
#		list(REMOVE_ITEM cls_sources ${f})
#	endif()
#endforeach()

# RPBA files
file(GLOB rpba_sources rpba.cpp rpbacore.cpp rbacore.cpp)

# SIM files
set(sim_sources sim.cpp)


set(stlplus_files stlplus/portability/file_system.cpp stlplus/portability/portability_fixes.cpp stlplus/portability/wildcard.cpp)
list(APPEND cls_sources ${stlplus_files})



set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "limited configs" FORCE)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build mode (Values: Release, Debug)")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)

# OpenCV
if(CMAKE_COMPILER_IS_GNUCXX)
	set(OPENCV_IMGCODECS_LIB "libopencv_imgcodecs.so")
	set(OPENCV_LIBS libopencv_core.so libopencv_imgproc.so libopencv_highgui.so libopencv_features2d.so)
else()
	set(OPENCV_VERSION_SUFFIX "300" CACHE STRING "OpenCV version suffix")
	set(opencv_libs opencv_core opencv_highgui opencv_imgproc)
	foreach(f ${opencv_libs})
		list(APPEND OPENCV_LIBS "${f}${OPENCV_VERSION_SUFFIX}.lib")
	endforeach()
	set(OPENCV_IMGCODECS_LIB "opencv_imgcodecs${OPENCV_VERSION_SUFFIX}.lib")
endif()
find_path(OPENCV_INCLUDE_DIR opencv/highgui.h opencv/cv.h $ENV{OPENCV_ROOT}/include)
find_path(OPENCV_LIB_DIR ${OPENCV_LIBS} $ENV{OPENCV_ROOT}/lib /usr/local/lib)
include_directories(${OPENCV_INCLUDE_DIR})
link_directories(${OPENCV_LIB_DIR})
set(CLS_TARGET_LIBS ${CLS_TARGET_LIBS} ${OPENCV_LIBS})
file(GLOB opencv_libs ${OPENCV_LIB_DIR}/*.*)
foreach(f ${opencv_libs})
	if( ${f} MATCHES ${OPENCV_IMGCODECS_LIB} )
		list(FIND OPENCV_LIBS ${OPENCV_IMGCODECS_LIB} index)
		if( ${index} LESS 0 )
			list(APPEND OPENCV_LIBS ${OPENCV_IMGCODECS_LIB})
			message(STATUS "OpenCV 3.x.x detected -> add module IMGCODECS")
		else()
			message(STATUS "OpenCV 3.x.x detected")
		endif()
	endif()
endforeach()
set(CLS_TARGET_LIBS ${CLS_TARGET_LIBS} ${OPENCV_LIBS})


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
			set(CLS_TARGET_LIBS ${CLS_TARGET_LIBS} tcmalloc_and_profiler)
			set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} tcmalloc_and_profiler)
		else()
			set(CLS_TARGET_LIBS ${CLS_TARGET_LIBS} tcmalloc)
			set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} tcmalloc)
		endif()
	endif()

	# TCMalloc with profiler
	if(BUILD_WITH_PROFILER)
		if(BUILD_WITH_TCMALLOC)
		else()
			set(CLS_TARGET_LIBS ${CLS_TARGET_LIBS} tcmalloc_and_profiler)
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
  set(CLS_TARGET_LIBS ${CLS_TARGET_LIBS} metis.so)
  set(rpba_TARGET_LIBS ${rpba_TARGET_LIBS} metis.so)
endif()
      
add_executable(cls ${cls_sources})
target_link_libraries(cls ${CLS_TARGET_LIBS})

add_executable(rpba ${rpba_sources})
target_link_libraries(rpba ${rpba_TARGET_LIBS})

add_executable(sim ${sim_sources})
target_link_libraries(sim ${SIM_TARGET_LIBS})
