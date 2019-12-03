
project ( waternetwork )

cmake_minimum_required( VERSION 2.8.11 )

# CGAL and its components
set( CGAL_DIR /localdisk/env/CGAL-5.0)
find_package( CGAL QUIET COMPONENTS )

if ( NOT CGAL_FOUND )

  message( STATUS "This project requires the CGAL library, and will not be compiled." )
  return()  

endif()

# Boost and its components
set ( BOOST_ROOT /localdisk/env/boost )
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

set( Boost_USE_MULTITHREADED ON )

################# include GROMACS
if ( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif()

# CMake modules are in a subdirectory to keep this file cleaner
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# In principle, this could be deduced from GROMACS_IS_DOUBLE returned by
# find_package(GROMACS) based on the suffix alone, but it is clearer that the
# user explicitly sets what they want to get, and then need to provide a suffix
# to match.
option( GMX_DOUBLE "Use double precision" OFF )
set( GMX_SUFFIX "_2019.4" CACHE STRING "Suffix for the GROMACS installation to use (empty for default)" )

# This does not allow for a non-suffixed double-precision libgromacs, but
# that should be rare enough for demonstration purposes.
if ( GMX_DOUBLE AND NOT GMX_SUFFIX )
    set( GROMACS_SUFFIX "_d" )
else()
    set(GROMACS_SUFFIX ${GMX_SUFFIX})
endif()

find_package(GROMACS 2019.4 REQUIRED)
gromacs_check_double(GMX_DOUBLE)
gromacs_check_compiler(CXX)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

set (CMAKE_CXX_FLAGS "-Wall -O2 -pg -pthread -std=c++17 -lfftw3 -fext-numeric-literals")

add_library( hbn SHARED
  src/cgal.cpp
  src/calculate-hydrogen-bond-energy.cpp
  src/brute-find-hydrogen-bonds.cpp
  src/delaunay-find-hydrogen-bonds.cpp
  src/alpha-shape-water-locator.cpp)
include_directories( hbn PUBLIC include )

include_directories( include )

set( SOURCE_WATERNETWORK src/waternetwork.cpp )

set( SOURCE_ALPHASHAPE src/alphashape.cpp )

set( SOURCE_TESTS test/tests.cpp )


# Creating entries for target: waternetwork
# ############################
add_executable( waternetwork ${SOURCE_WATERNETWORK} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS waternetwork )

# Link the executable to CGAL and third-party libraries
target_link_libraries(waternetwork ${GROMACS_LIBRARIES} CGAL::CGAL hbn)


# Creating entries for target: alphashape
# ############################
add_executable( alphashape ${SOURCE_ALPHASHAPE} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS alphashape )

# Link the executable to CGAL and third-party libraries
target_link_libraries(alphashape ${GROMACS_LIBRARIES} CGAL::CGAL hbn)


# Creating entries for target: tests
# ############################
add_executable( tests ${SOURCE_TESTS} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS tests )

# Link the executable to CGAL and third-party libraries
target_link_libraries(tests ${GROMACS_LIBRARIES} CGAL::CGAL hbn fftw3 )

enable_testing()
add_test(NAME tests COMMAND tests)

