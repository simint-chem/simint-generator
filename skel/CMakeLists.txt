cmake_minimum_required(VERSION 3.0.2)
project(simint C)
enable_language(CXX)

##################################################
# Flags passed in
##################################################
set(SIMINT_VECTOR  "" CACHE BOOL "Type of vectorization to use")


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


# Install the config file
configure_file("simintConfig.cmake.in" "simintConfig.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/simintConfig.cmake DESTINATION share/cmake/simint)

# Configure & install the vector_config header
configure_file("simint/vectorization/vector_config.h.in"
               "${CMAKE_BINARY_DIR}/stage/simint/vectorization/vector_config.h"
               @ONLY)

foreach(def ${SIMINT_CONFIG_DEFINES})
  file(APPEND "${CMAKE_BINARY_DIR}/stage/simint/vectorization/vector_config.h"
              "#define ${def}\n")
endforeach()

install(FILES ${CMAKE_BINARY_DIR}/stage/simint/vectorization/vector_config.h DESTINATION include/simint/vectorization)


# Main code subdirectory
add_subdirectory(simint)

# Make the test executables
add_subdirectory(test)

# Build examples
add_subdirectory(examples)
