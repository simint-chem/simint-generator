cmake_minimum_required(VERSION 3.0)
project(simint C)
enable_language(Fortran) # for liberd

set(SIMINT_SIMD FALSE CACHE BOOL "Enable SIMD")

list(APPEND SIMINT_C_INCLUDES "${CMAKE_SOURCE_DIR}")

list(APPEND SIMINT_C_FLAGS "-std=c99")
list(APPEND SIMINT_C_FLAGS "-Wall -Wunused-variable -pedantic")
list(APPEND SIMINT_C_FLAGS "-Werror -pedantic-errors")
list(APPEND SIMINT_C_FLAGS "-ggdb")

list(APPEND SIMINT_LINK_FLAGS "-Wall")   #??


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  list(APPEND SIMINT_C_FLAGS "-wd10397 -xavx -qopt-report=4")
endif()


# Check for OpenMP for pragma omp simd lines
find_package(OpenMP REQUIRED)
list(APPEND SIMINT_C_FLAGS "-DSIMINT_SIMD")
list(APPEND SIMINT_C_FLAGS "${OpenMP_C_FLAGS}")
list(APPEND SIMINT_LINK_FLAGS "${OpenMP_C_FLAGS}")


# Turn lists into a string
# Except for includes. They are handled differently in cmake...
string(REPLACE ";" " " SIMINT_C_FLAGS "${SIMINT_C_FLAGS}")
string(REPLACE ";" " " SIMINT_LINK_FLAGS "${SIMINT_LINK_FLAGS}")

message(STATUS "   C_FLAGS: ${SIMINT_C_FLAGS}")
message(STATUS "C_INCLUDES: ${SIMINT_C_INCLUDES}")
message(STATUS "LINK_FLAGS: ${SIMINT_LINK_FLAGS}")


# Add the subdirectories with the libraries test executables
add_subdirectory(boys)
add_subdirectory(eri)
add_subdirectory(valeev)
add_subdirectory(test)





##################################################
# LibERD
##################################################
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  list(APPEND LIBERD_F90_FLAGS "-Wall -g")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  list(APPEND LIBERD_F90_FLAGS "-warn all -g")
else()  
  list(APPEND LIBERD_F90_FLAGS "-Wall -g")
endif() 

list(APPEND LIBERD_F90_FLAGS "${OpenMP_CXX_FLAGS}")
add_subdirectory(liberd)
