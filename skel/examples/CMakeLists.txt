###########################
# Add examples
###########################
add_executable(example1 example1.c)
target_compile_options(example1 PRIVATE ${SIMINT_C_FLAGS})
target_include_directories(example1 PRIVATE ${SIMINT_C_INCLUDES})
target_link_libraries(example1 PRIVATE simint m ${SIMINT_C_LINK_FLAGS})

add_executable(example_packed example_packed.c)
target_compile_options(example_packed PRIVATE ${SIMINT_C_FLAGS})
target_include_directories(example_packed PRIVATE ${SIMINT_C_INCLUDES})
target_link_libraries(example_packed PRIVATE simint m ${SIMINT_C_LINK_FLAGS})

