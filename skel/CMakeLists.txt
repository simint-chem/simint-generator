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


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5 -wd10397 -wd981")
  #list(APPEND SIMINT_C_FLAGS "-Werror -pedantic-errors")
else()
  message(FATAL_ERROR "Unsupported compiler (for now)")
endif()


##################################################
# AVX, etc
##################################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if(SIMINT_SSE)
    list(APPEND SIMINT_C_FLAGS "-msse2 -DSIMINT_SIMD_LEN=2")
  elseif(SIMINT_AVX)
    list(APPEND SIMINT_C_FLAGS "-xavx -DSIMINT_SIMD_LEN=4")
  elseif(SIMINT_MIC)
    list(APPEND SIMINT_C_FLAGS "-mmic -DSIMINT_SIMD_LEN=8")
  endif()

  # Separate, since we may want -mmic and -no-vec, etc
  if(SIMINT_SCALAR)
    list(APPEND SIMINT_C_FLAGS "-no-vec -no-simd -DSIMINT_SIMD_LEN=1")
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

