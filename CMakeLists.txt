# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 2.8)

project ( waternetwork )

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
set(GROMACS_DIR "/localdisk/local/share/cmake/gromacs")

# This does not allow for a non-suffixed double-precision libgromacs, but
# that should be rare enough for demonstration purposes.
if (GMX_DOUBLE AND NOT GMX_SUFFIX)
    set(GROMACS_SUFFIX "_d")
else()
    set(GROMACS_SUFFIX ${GMX_SUFFIX})
endif()

find_package(GROMACS 5.1 REQUIRED)
gromacs_check_double(GMX_DOUBLE)
gromacs_check_compiler(CXX)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

set (CMAKE_CXX_FLAGS "-Wall -O2 -pg -std=c++14")

# includes and sources
include_directories(include/ include/waternetwork include/alphashape)
#file(GLOB SOURCES "src/*.cpp")

set( SOURCE_WATERNETWORK src/waternetwork.cpp
  src/dipole_module.cpp
  src/graph_module.cpp
  src/lifetime_module.cpp)

set( SOURCE_ALPHASHAPE src/alphashape.cpp
  src/alpha_shape_module.cpp)

# Creating entries for all C++ files with "main" routine
# ##########################################################

add_executable ( waternetwork ${SOURCE_WATERNETWORK})
add_executable ( alphashape ${SOURCE_ALPHASHAPE})

target_link_libraries( waternetwork ${GROMACS_LIBRARIES} )
target_link_libraries( alphashape ${GROMACS_LIBRARIES}
  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )



