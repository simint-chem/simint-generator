cmake_minimum_required(VERSION 3.0)
project(simint C)
enable_language(Fortran) # for liberd

set(SIMINT_AVX FALSE CACHE BOOL "Enable AVX")
set(SIMINT_MIC FALSE CACHE BOOL "Enable Xeon Phi")

list(APPEND SIMINT_C_INCLUDES "${CMAKE_SOURCE_DIR}")


##################################################
# Basic options for simint
##################################################
list(APPEND SIMINT_C_FLAGS "-std=c99")
list(APPEND SIMINT_C_FLAGS "-O3")
list(APPEND SIMINT_C_FLAGS "-Wall -Wunused-variable -pedantic")
list(APPEND SIMINT_C_FLAGS "-Werror -pedantic-errors")
list(APPEND SIMINT_C_FLAGS "-ggdb")

list(APPEND SIMINT_LINK_FLAGS "-Wall")   #??


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=4 -wd10397")
endif()


##################################################
# Basic options for libERD
##################################################
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  list(APPEND LIBERD_F90_FLAGS "-Wall -g")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  list(APPEND LIBERD_F90_FLAGS "-warn all -g")
else()  
  list(APPEND LIBERD_F90_FLAGS "-Wall -g")
endif() 




##################################################
# AVX, etc
##################################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if(SIMINT_AVX)
    list(APPEND SIMINT_C_FLAGS "-xavx")
    list(APPEND SIMINT_LINK_FLAGS "-xavx")
  elseif(SIMINT_MIC)
    list(APPEND SIMINT_C_FLAGS "-mmic")
    list(APPEND SIMINT_LINK_FLAGS "-mmic")
  endif()
endif()

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  if(SIMINT_AVX)
    list(APPEND LIBERD_F90_FLAGS "-xavx")
  elseif(SIMINT_MIC)
    list(APPEND LIBERD_F90_FLAGS "-mmic")
  endif()
endif()


# Turn lists into a string
# Except for includes. They are handled differently in cmake...
string(REPLACE ";" " " SIMINT_C_FLAGS "${SIMINT_C_FLAGS}")
string(REPLACE ";" " " SIMINT_LINK_FLAGS "${SIMINT_LINK_FLAGS}")
string(REPLACE ";" " " LIBERD_F90_FLAGS "${LIBERD_F90_FLAGS}")

message(STATUS "  C_INCLUDES: ${SIMINT_C_INCLUDES}")
message(STATUS "     Simint C_FLAGS: ${SIMINT_C_FLAGS}")
message(STATUS "Simint C LINK_FLAGS: ${SIMINT_LINK_FLAGS}")
message(STATUS "   liberd F90 FLAGS: ${LIBERD_F90_FLAGS}")


# Add the subdirectories with the libraries test executables
add_subdirectory(liberd)
add_subdirectory(boys)
add_subdirectory(eri)
add_subdirectory(valeev)
add_subdirectory(test)



