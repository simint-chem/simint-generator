###################
# Generator flags
###################
list(APPEND GENERATOR_CXX_FLAGS "-std=c++11")
list(APPEND GENERATOR_CXX_FLAGS "-O3")
list(APPEND GENERATOR_CXX_FLAGS "-Wall -Wunused-variable -pedantic")
#list(APPEND GENERATOR_CXX_FLAGS "-Werror -pedantic-errors")
list(APPEND GENERATOR_CXX_FLAGS "-ggdb")

list(APPEND GENERATOR_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")


add_executable(generator generator.cpp
                         AlgorithmBase.cpp
                         Helpers.cpp
                         Boys.cpp
                         WriterBase.cpp
                         HRRWriter.cpp
                         ETWriter.cpp
                         VRRWriter.cpp
                         FileWriter.cpp
              )

string(REPLACE ";" " " GENERATOR_CXX_FLAGS "${GENERATOR_CXX_FLAGS}")
set_target_properties(generator PROPERTIES INCLUDE_DIRECTORIES "${GENERATOR_CXX_INCLUDES}")
set_target_properties(generator PROPERTIES COMPILE_FLAGS ${GENERATOR_CXX_FLAGS})

# No special linking needed
