###############################
# Builds the Simint Generator
###############################


###################
# Options
###################
set(SIMINT_ORDERING "Psi4" CACHE STRING "Ordering of the resulting integrals")
set_property(CACHE SIMINT_ORDERING PROPERTY STRINGS Psi4 GAMESS) 

###################
# Generator flags
###################
list(APPEND GENERATOR_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")


if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   383 : value copied to temporary, reference to temporary used
  #   869 : parameter <var> never referenced
  #         lots of false positives due to intel bug
  #  1418 : function definition with no prior declaration
  # 11074 : Inlining inhibited
  # 11076 : Use qopt-report for inline stuff
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-w3;-wd981;-wd383;-wd869;-wd1418;-wd11074;-wd11076")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "PGI")
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-Wall")
else()
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-Wall;-Wextra;-pedantic;-Wshadow;-Wno-unused-parameter")
endif()


#################################
# Determine the output ordering
#################################
if("${SIMINT_ORDERING}" MATCHES "Psi4")
  set(SIMINT_ORDERING_FILE Ordering_Psi4.cpp)
elseif("${SIMINT_ORDERING}" MATCHES "GAMESS")
  set(SIMINT_ORDERING_FILE Ordering_GAMESS.cpp)
else()
  message(FATAL_ERROR "Unknown ordering ${SIMINT_ORDERING}")
endif()



#################################
# Common generator functionality
#################################
add_library(generator_common OBJECT
                         CommandLine.cpp
                         Types.cpp

                         ${SIMINT_ORDERING_FILE}
           )


target_include_directories(generator_common PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(generator_common PRIVATE ${GENERATOR_CXX_FLAGS})


#####################################
# Obara-Saika Two-electron Integrals
#####################################
add_subdirectory(ostei)
