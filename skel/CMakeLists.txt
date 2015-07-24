cmake_minimum_required(VERSION 2.8.10)
project(simint C)
enable_language(CXX)

###################
# OPTIONS
###################
set(SIMINT_AVX FALSE CACHE BOOL "Enable AVX")
set(SIMINT_MIC FALSE CACHE BOOL "Enable Xeon Phi")


##################################################
# Basic options for simint
##################################################
list(APPEND SIMINT_C_INCLUDES "${CMAKE_SOURCE_DIR}")
list(APPEND SIMINT_C_FLAGS "-std=c99")
list(APPEND SIMINT_C_FLAGS "-O3")
list(APPEND SIMINT_C_FLAGS "-Wall -Wunused-variable -pedantic")
#list(APPEND SIMINT_C_FLAGS "-Werror -pedantic-errors")
list(APPEND SIMINT_C_FLAGS "-ggdb")


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5 -wd10397 -vec-threshold0 -override-limits")
endif()


##################################################
# AVX, etc
##################################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if(SIMINT_AVX)
    list(APPEND SIMINT_C_FLAGS "-xavx")
  elseif(SIMINT_SSE)
    list(APPEND SIMINT_C_FLAGS "-msse2")
  elseif(SIMINT_MIC)
    list(APPEND SIMINT_C_FLAGS "-mmic")
  endif()

  # Separate, since we may want -mmic and -no-vec, etc
  if(SIMINT_SCALAR)
    list(APPEND SIMINT_C_FLAGS "-no-vec -no-simd")
  endif()

endif()


##################################################
# OpenMP dependency
##################################################
find_package(OpenMP REQUIRED)
list(APPEND SIMINT_C_FLAGS "${OpenMP_C_FLAGS}") 


# Turn lists into a string
# Except for includes. They are handled differently in cmake...
string(REPLACE ";" " " SIMINT_C_FLAGS "${SIMINT_C_FLAGS}")

message(STATUS "  C_INCLUDES: ${SIMINT_C_INCLUDES}")
message(STATUS "     Simint C_FLAGS: ${SIMINT_C_FLAGS}")


# Add the subdirectories with the libraries test executables
add_subdirectory(boys)
add_subdirectory(shell)
add_subdirectory(eri)
add_subdirectory(test)
