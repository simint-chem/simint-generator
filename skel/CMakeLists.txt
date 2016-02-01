cmake_minimum_required(VERSION 3.0.2)
project(simint C)
enable_language(CXX)


##################################################
# Basic options for simint
##################################################
list(APPEND SIMINT_C_INCLUDES "${CMAKE_SOURCE_DIR}")


##################################################
# Flags passed in
##################################################
set(SIMINT_C_FLAGS CACHE BOOL "Additional compilation flags for Simint")
set(SIMINT_CXX_FLAGS CACHE BOOL "Additional compilation flags for Simint testing/benchmarking")



if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5;-w3;-wd10397;-wd981")
else()
  message(FATAL_ERROR "Unsupported compiler (for now)")
endif()


##################################################
# OpenMP dependency
##################################################
find_package(OpenMP REQUIRED)
list(APPEND SIMINT_C_FLAGS "${OpenMP_C_FLAGS}") 


message(STATUS "  C_INCLUDES: ${SIMINT_C_INCLUDES}")
message(STATUS "     Simint C_FLAGS: ${SIMINT_C_FLAGS}")



# Add the subdirectories with the library components
add_subdirectory(boys)
add_subdirectory(shell)
add_subdirectory(eri)



# Create the main library
add_library(simint STATIC
                   simint_init.c
                   $<TARGET_OBJECTS:simint_eri>
                   $<TARGET_OBJECTS:simint_shell>
                   $<TARGET_OBJECTS:simint_boys>
           )


# Make the test executables
add_subdirectory(test)
