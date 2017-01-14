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
list(APPEND SIMINT_C_INCLUDES "${CMAKE_BINARY_DIR}/stage") # some CMake-generated files go here

# Default to "Release"
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release")
endif()

##################################################
# Read the default flags and generated flags file
##################################################
include(cmake/DefaultFlags.cmake)


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

# Install the config file
configure_file("simintConfig.cmake.in" "simintConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/simintConfig.cmake DESTINATION share/cmake/simint)

# Configure & install the simint_config header
configure_file("simint/simint_config.h.in" "${CMAKE_BINARY_DIR}/stage/simint/simint_config.h" @ONLY)

foreach(def ${SIMINT_CONFIG_DEFINES})
  file(APPEND "${CMAKE_BINARY_DIR}/stage/simint/simint_config.h" "#define ${def}\n")
endforeach()

install(FILES ${CMAKE_BINARY_DIR}/stage/simint/simint_config.h DESTINATION include/simint)

# Main code subdirectory
add_subdirectory(simint)

# Make the test executables
add_subdirectory(test)

# Build examples
add_subdirectory(examples)
