#######################
# Common OS OSTEI stuff
#######################
add_library(ostei_generator_common OBJECT
                         OSTEI_VRR_Algorithm_Base.cpp
                         OSTEI_HRR_Algorithm_Base.cpp
                         OSTEI_VRR_Writer.cpp
                         OSTEI_HRR_Writer.cpp
                         OSTEI_Writer.cpp
           )


target_include_directories(ostei_generator_common PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(ostei_generator_common PRIVATE ${GENERATOR_CXX_FLAGS})


add_executable(ostei_generator ostei_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:ostei_generator_common>)
target_include_directories(ostei_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(ostei_generator PRIVATE ${GENERATOR_CXX_FLAGS})

add_executable(ostei_vrr_generator ostei_vrr_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:ostei_generator_common>)
target_include_directories(ostei_vrr_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(ostei_vrr_generator PRIVATE ${GENERATOR_CXX_FLAGS})

add_executable(ostei_hrr_generator ostei_hrr_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:ostei_generator_common>)
target_include_directories(ostei_hrr_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(ostei_hrr_generator PRIVATE ${GENERATOR_CXX_FLAGS})
