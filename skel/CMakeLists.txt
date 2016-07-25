cmake_minimum_required(VERSION 3.0.2)
project(simint C)
enable_language(CXX)

##################################################
# Flags passed in
##################################################
set(SIMINT_C_FLAGS CACHE BOOL "Additional compilation flags for Simint")
set(SIMINT_CXX_FLAGS CACHE BOOL "Additional compilation flags for Simint testing/benchmarking")
set(SIMINT_BUILD_SHARED FALSE CACHE BOOL "Build a shared library, rather than a static library")


##################################################
# Basic options for simint
##################################################
list(APPEND SIMINT_C_INCLUDES "${CMAKE_SOURCE_DIR}")


##################################################
# Flags for warnings (and disabling warnings)
##################################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   869 : parameter "xxx" never referenced (ok, but try compiling without it sometimes)
  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_C_FLAGS "-wd10397;-wd981;-wd869")

  list(APPEND SIMINT_CXX_FLAGS "-std=c++11;-restrict")
  list(APPEND SIMINT_CXX_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_CXX_FLAGS "-wd10397;-wd981;-wd869")
elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-parameter")

  list(APPEND SIMINT_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_CXX_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_CXX_FLAGS "-Wno-unused-parameter")
else()
  message(FATAL_ERROR "Unsupported compiler")
endif()


##################################################
# OpenMP dependency
##################################################
find_package(OpenMP REQUIRED)
list(APPEND SIMINT_C_FLAGS "${OpenMP_C_FLAGS}") 
list(APPEND SIMINT_CXX_FLAGS "${OpenMP_C_FLAGS}") 


message(STATUS "  Simint C_INCLUDES: ${SIMINT_C_INCLUDES}")
message(STATUS "     Simint C_FLAGS: ${SIMINT_C_FLAGS}")
message(STATUS "   Simint CXX_FLAGS: ${SIMINT_C_FLAGS}")


if(${SIMINT_BUILD_SHARED})
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
