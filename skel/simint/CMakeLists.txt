# This cmake files builds and installs the main simint library
# and its header files


# Find all the eri files
file(GLOB SIMINT_ERI_FILES eri/gen/*.c)

# All the files needed for the build
set(SIMINT_LIBRARY_FILES ${SIMINT_ERI_FILES}

                         simint_init.c

                         boys/boys_longfac.c
                         boys/boys_shortgrid.c

                         eri/eri_init.c
                         eri/eri.c

                         shell/shell.c
                         shell/shell_constants.c
    )


# Create the main library
add_library(simint ${SIMINT_LIBRARY_TYPE} ${SIMINT_LIBRARY_FILES})
target_include_directories(simint PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint PRIVATE $<$<COMPILE_LANGUAGE:C>:${SIMINT_C_FLAGS}>)

install(TARGETS simint ARCHIVE DESTINATION lib
                       LIBRARY DESTINATION lib)

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" DESTINATION include FILES_MATCHING PATTERN "*.h")
