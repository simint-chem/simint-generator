###################
# Generator flags
###################
list(APPEND GENERATOR_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   383 : value copied to temporary, reference to temporary used
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-qopt-report=5;-w3;-wd10397;-wd981;-wd383")
else()
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-Wall;-Wextra;-pedantic;-Wno-unused-parameter")
endif()



add_library(generator_common OBJECT
                         BoysGenerator.cpp
                         BoysGenerator_FO.cpp
                         BoysGenerator_split.cpp

                         GeneratorInfoBase.cpp
                         Types.cpp
                         Naming.cpp
                         Printing.cpp
                         CommandLine.cpp
           )


target_include_directories(generator_common PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(generator_common PRIVATE ${GENERATOR_CXX_FLAGS})
set_target_properties(generator_common PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(generator_common PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})


add_subdirectory(eri)
