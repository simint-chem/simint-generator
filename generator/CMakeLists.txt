###################
# Generator flags
###################
list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
list(APPEND GENERATOR_CXX_FLAGS "-O3")
list(APPEND GENERATOR_CXX_FLAGS "-Wall -Wunused-variable -pedantic")
#list(APPEND GENERATOR_CXX_FLAGS "-Werror -pedantic-errors")
list(APPEND GENERATOR_CXX_FLAGS "-ggdb")

list(APPEND GENERATOR_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")

string(REPLACE ";" " " GENERATOR_CXX_FLAGS "${GENERATOR_CXX_FLAGS}")

add_library(generator_common STATIC
                         VRR_Algorithm_Base.cpp
                         ET_Algorithm_Base.cpp
                         HRR_Algorithm_Base.cpp
                         Helpers.cpp
                         Boys.cpp
                         WriterInfo.cpp
                         VRR_Writer.cpp
                         ET_Writer.cpp
                         HRR_Writer.cpp
                         FileWriter.cpp
           )
set_target_properties(generator_common PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(generator_common PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})


add_executable(eri_generator eri_generator.cpp)
set_target_properties(eri_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(eri_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})
target_link_libraries(eri_generator generator_common)


add_executable(vrr_generator vrr_generator.cpp)
set_target_properties(vrr_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(vrr_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})
target_link_libraries(vrr_generator generator_common)

add_executable(hrr_generator hrr_generator.cpp)
set_target_properties(hrr_generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(hrr_generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})
target_link_libraries(hrr_generator generator_common)
