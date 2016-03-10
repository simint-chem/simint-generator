###############################
# Builds the Simint Generator
###############################


###################
# Generator flags
###################
list(APPEND GENERATOR_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
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
else()
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-Wall;-Wextra;-pedantic;-Wshadow;-Wno-unused-parameter")
endif()


#################################
# Common generator functionality
#################################
add_library(generator_common OBJECT
                         BoysGenerator.cpp
                         BoysGenerator_FO.cpp
                         BoysGenerator_split.cpp

                         GeneratorInfoBase.cpp
                         Types.cpp
                         CommandLine.cpp
           )


target_include_directories(generator_common PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(generator_common PRIVATE ${GENERATOR_CXX_FLAGS})


#######################
# Obara-Saika ERI
#######################
add_subdirectory(eri)
