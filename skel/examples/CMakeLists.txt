###########################
# Add examples
###########################
add_executable(example1 example1.c)
target_compile_options(example1 PRIVATE ${SIMINT_C_FLAGS})
target_include_directories(example1 PRIVATE ${SIMINT_C_INCLUDES})
target_link_libraries(example1 PRIVATE simint m ${SIMINT_LINK_FLAGS})

add_executable(example_packed example_packed.c)
target_compile_options(example_packed PRIVATE ${SIMINT_C_FLAGS})
target_include_directories(example_packed PRIVATE ${SIMINT_C_INCLUDES})
target_link_libraries(example_packed PRIVATE simint m ${SIMINT_LINK_FLAGS})

add_executable(example_water_spsp example_water_spsp.c)
target_compile_options(example_water_spsp PRIVATE ${SIMINT_C_FLAGS})
target_include_directories(example_water_spsp PRIVATE ${SIMINT_C_INCLUDES})
target_link_libraries(example_water_spsp PRIVATE simint m ${SIMINT_LINK_FLAGS})

if(ENABLE_FORTRAN)
    add_executable(example_water_spsp_f90 example_water_spsp.f90)
    target_compile_options(example_water_spsp_f90 PRIVATE ${SIMINT_Fortran_FLAGS})
    target_include_directories(example_water_spsp_f90 PRIVATE ${SIMINT_Fortran_INCLUDES})
    target_link_libraries(example_water_spsp_f90 PRIVATE simint m ${SIMINT_LINK_FLAGS})
endif()

if(${simint_MAXDER} GREATER 0)
    add_executable(example_water_ppps_deriv example_water_ppps_deriv.c)
    target_compile_options(example_water_ppps_deriv PRIVATE ${SIMINT_C_FLAGS})
    target_include_directories(example_water_ppps_deriv PRIVATE ${SIMINT_C_INCLUDES})
    target_link_libraries(example_water_ppps_deriv PRIVATE simint m ${SIMINT_LINK_FLAGS})
endif()

