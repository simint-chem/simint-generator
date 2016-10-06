cmake_minimum_required(VERSION 3.0.2)
project(simint C)
enable_language(CXX)

##################################################
# Flags passed in
##################################################
set(SIMINT_C_FLAGS         CACHE STRING "Additional compilation flags for Simint")
set(SIMINT_C_LINK_FLAGS    CACHE STRING "Additional linker flags for the Simint library")
set(SIMINT_TESTS_CXX_FLAGS CACHE STRING "Additional compilation flags for Simint testing/benchmarking")
set(SIMINT_BUILD_SHARED    FALSE CACHE BOOL "Build a shared library, rather than a static library")


##################################################
# Basic options for simint
##################################################
list(APPEND SIMINT_C_INCLUDES "${CMAKE_SOURCE_DIR}")

# Default to "Release"
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release")
endif()

##################################################
# Read the default flags and generated flags file
##################################################
include(cmake/DefaultFlags.cmake)
include(cmake/DefaultFlags_generated.cmake)


##################################################
# OpenMP dependency
##################################################
find_package(OpenMP REQUIRED)
list(APPEND SIMINT_TESTS_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 


message(STATUS "    SIMINT_C_INCLUDES: ${SIMINT_C_INCLUDES}")
message(STATUS "       SIMINT_C_FLAGS: ${SIMINT_C_FLAGS}")
message(STATUS "  SIMINT_C_LINK_FLAGS: ${SIMINT_C_LINK_FLAGS}")


if(SIMINT_BUILD_SHARED)
  set(SIMINT_LIBRARY_TYPE "SHARED")
  message(STATUS "Building shared library")
else()
  set(SIMINT_LIBRARY_TYPE "STATIC")
  message(STATUS "Building static library")
endif()



# Main code subdirectory
add_subdirectory(simint)

# Make the test executables
add_subdirectory(test)

# Build examples
add_subdirectory(examples)
