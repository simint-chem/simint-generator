# Add the subdirectories with the library components
add_subdirectory(boys)
add_subdirectory(shell)
add_subdirectory(eri)


# Create the main library
add_library(simint STATIC
                   simint_init.c
                   $<TARGET_OBJECTS:simint_eri>
                   $<TARGET_OBJECTS:simint_shell>
                   $<TARGET_OBJECTS:simint_boys>
           )
target_include_directories(simint PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint PRIVATE ${SIMINT_C_FLAGS})

install(TARGETS simint ARCHIVE DESTINATION lib)


