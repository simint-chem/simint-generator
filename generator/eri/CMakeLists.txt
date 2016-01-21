add_library(eri_generator_common OBJECT
                         ERI_VRR_Algorithm_Base.cpp
                         ERI_ET_Algorithm_Base.cpp
                         ERI_HRR_Algorithm_Base.cpp
                         ERI_VRR_Writer.cpp
                         ERI_ET_Writer.cpp
                         ERI_HRR_Writer.cpp
                         ERI_Writer.cpp
           )


target_include_directories(eri_generator_common PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(eri_generator_common PRIVATE ${GENERATOR_CXX_FLAGS})


add_executable(eri_generator eri_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:eri_generator_common>)
target_include_directories(eri_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(eri_generator PRIVATE ${GENERATOR_CXX_FLAGS})

add_executable(eri_vrr_generator eri_vrr_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:eri_generator_common>)
target_include_directories(eri_vrr_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(eri_vrr_generator PRIVATE ${GENERATOR_CXX_FLAGS})

add_executable(eri_et_generator eri_et_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:eri_generator_common>)
target_include_directories(eri_et_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(eri_et_generator PRIVATE ${GENERATOR_CXX_FLAGS})

add_executable(eri_hrr_generator eri_hrr_generator.cpp $<TARGET_OBJECTS:generator_common> $<TARGET_OBJECTS:eri_generator_common>)
target_include_directories(eri_hrr_generator PRIVATE ${GENERATOR_CXX_INCLUDES})
target_compile_options(eri_hrr_generator PRIVATE ${GENERATOR_CXX_FLAGS})
