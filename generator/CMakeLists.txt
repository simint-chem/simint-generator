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
  list(APPEND GENERATOR_CXX_FLAGS "-qopt-report=5 -w3 -wd10397 -wd981 -wd383")
  list(APPEND GENERATOR_CXX_FLAGS "-O0 -ggdb")
else()
  list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
  list(APPEND GENERATOR_CXX_FLAGS "-Wall -Wextra -pedantic -Wno-unused-parameter")
  list(APPEND GENERATOR_CXX_FLAGS "-O0 -ggdb")
endif()



# Prepare list for CXX Flags
string(REPLACE ";" " " GENERATOR_CXX_FLAGS "${GENERATOR_CXX_FLAGS}")

add_library(generator_common OBJECT
                         Classes.cpp
                         Naming.cpp
                         Printing.cpp
                         CommandLine.cpp
                         BoysGenerator.cpp
                         BoysGenerator_FO.cpp
                         BoysGenerator_split.cpp

                         ERIGeneratorInfo.cpp
                         ERI_VRR_Algorithm_Base.cpp
                         ERI_ET_Algorithm_Base.cpp
                         ERI_HRR_Algorithm_Base.cpp
                         ERI_VRR_Writer.cpp
                         ERI_ET_Writer.cpp
                         ERI_HRR_Writer.cpp
                         ERI_Writer.cpp
           )



set_target_properties(generator_common PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(generator_common PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})



add_executable(eri_generator eri_generator.cpp $<TARGET_OBJECTS:generator_common>)
set_target_properties(eri_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(eri_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})

add_executable(vrr_generator vrr_generator.cpp $<TARGET_OBJECTS:generator_common>)
set_target_properties(vrr_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(vrr_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})

add_executable(et_generator et_generator.cpp $<TARGET_OBJECTS:generator_common>)
set_target_properties(et_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(et_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})

add_executable(hrr_generator hrr_generator.cpp $<TARGET_OBJECTS:generator_common>)
set_target_properties(hrr_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(hrr_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})
