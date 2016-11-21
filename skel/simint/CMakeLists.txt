# This cmake files builds and installs the main simint library
# and its header files


# Find all the generated files
file(GLOB SIMINT_OSTEI_FILES ostei/gen/*.c)

# All the files needed for the build
set(SIMINT_LIBRARY_FILES ${SIMINT_OSTEI_FILES}

                         simint_init.c
                         simint_eri.c

                         boys/boys_shortgrid.c

                         ostei/ostei_init.c

                         ostei/ostei_general_vrr.c
                         ostei/ostei_general_hrr.c
                         ostei/ostei.c
                         ostei/recur_lookup.c

                         shell/shell.c
                         shell/shell_screen.c
                         shell/shell_constants.c
    )


# Create the main library
add_library(simint ${SIMINT_LIBRARY_TYPE} ${SIMINT_LIBRARY_FILES})
target_include_directories(simint PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint PRIVATE ${SIMINT_C_FLAGS})
target_link_libraries(simint PRIVATE ${SIMINT_C_LINK_FLAGS})

install(TARGETS simint ARCHIVE DESTINATION lib
                       LIBRARY DESTINATION lib)

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" DESTINATION include FILES_MATCHING PATTERN "*.h")
