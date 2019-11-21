# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 2.8)

project ( waternetwork 
  DESCRIPTION "Tools for analysis of hydrogen bonds networks from Molecular Dynamics Trajectories"
  LANGUAGES CXX )

set(Boost_USE_MULTITHREADED ON)
# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )

# include GROMACS
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# CMake modules are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# In principle, this could be deduced from GROMACS_IS_DOUBLE returned by
# find_package(GROMACS) based on the suffix alone, but it is clearer that the
# user explicitly sets what they want to get, and then need to provide a suffix
# to match.
option(GMX_DOUBLE "Use double precision" OFF)
set(GMX_SUFFIX "" CACHE STRING "Suffix for the GROMACS installation to use (empty for default)")
set(GROMACS_DIR "/usr/local/gromacs/2016.5")
#set(GROMACS_DIR "/localdisk/local/share/cmake/gromacs")

# This does not allow for a non-suffixed double-precision libgromacs, but
# that should be rare enough for demonstration purposes.
if (GMX_DOUBLE AND NOT GMX_SUFFIX)
    set(GROMACS_SUFFIX "_d")
else()
    set(GROMACS_SUFFIX ${GMX_SUFFIX})
endif()

find_package(GROMACS 2016.5 REQUIRED)
gromacs_check_double(GMX_DOUBLE)
gromacs_check_compiler(CXX)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

set (CMAKE_CXX_FLAGS "-fopenmp -Wall -O2 -pg -pthread -std=c++17 -lfftw3")

include_directories( include )

add_library ( cgal SHARED src/Cgal.cpp )
include_directories( cgal PUBLIC include )

set( SOURCE_WATERNETWORK src/waternetwork.cpp )

set( SOURCE_ALPHASHAPE src/alphashape.cpp )

set( SOURCE_TESTS test/tests.cpp )


add_executable ( waternetwork ${SOURCE_WATERNETWORK})
target_link_libraries( waternetwork ${GROMACS_LIBRARIES} cgal )

add_executable ( alphashape ${SOURCE_ALPHASHAPE})
target_link_libraries( alphashape ${GROMACS_LIBRARIES} cgal )

add_executable ( tests ${SOURCE_TESTS} )
target_link_libraries( tests ${GROMACS_LIBRARIES} cgal fftw3)

enable_testing()
add_test(NAME tests COMMAND tests)

