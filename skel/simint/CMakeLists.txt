# This cmake files builds and installs the main simint library
# and its header files


# Find which files we need to build
if(${SIMINT_MAXAM} GREATER ${SIMINT_MAXAM_AVAILABLE})
  message(FATAL_ERROR "Maximum angular momentum set is greater than I have available (${SIMINT_MAXAM_AVAILABLE})")
endif()

if(${SIMINT_MAXDER} GREATER ${SIMINT_MAXDER_AVAILABLE})
  message(FATAL_ERROR "Maximum derivative set is greater than I have available (${SIMINT_MAXDER_AVAILABLE})")
endif()


# Read the file list and convert to a cmake list
file(READ "ostei/gen/filelist_d${SIMINT_MAXDER}_${SIMINT_MAXAM}" SIMINT_OSTEI_FILES)
string(REGEX REPLACE "\n" ";" SIMINT_OSTEI_FILES "${SIMINT_OSTEI_FILES}")

# Configure the ostei_config.h file
configure_file("ostei/ostei_config.h.in"
               "${CMAKE_BINARY_DIR}/stage/simint/ostei/ostei_config.h"
               @ONLY)

# All the files needed for the build
set(SIMINT_LIBRARY_FILES ${SIMINT_OSTEI_FILES}

                         simint_init.c
                         simint_eri.c

                         boys/boys_shortgrid.c

                         ostei/ostei_init.c
                         ostei/ostei_deriv1_init.c

                         ostei/ostei_general_vrr.c
                         ostei/ostei_general_hrr.c
                         ostei/ostei.c
                         ostei/recur_lookup.c

                         shell/shell.c
                         shell/shell_screen.c
                         shell/shell_constants.c
    )

if(ENABLE_FORTRAN)
    set(SIMINT_LIBRARY_FILES ${SIMINT_LIBRARY_FILES} simint.f90)
    list(APPEND SIMINT_Fortran_FLAGS "-DSIMINT_OSTEI_MADXER=${SIMINT_MAXDER}")
endif()


# Create the main library
add_library(simint ${SIMINT_LIBRARY_TYPE} ${SIMINT_LIBRARY_FILES})
target_include_directories(simint PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint PRIVATE $<$<COMPILE_LANGUAGE:C>:${SIMINT_C_FLAGS}>
                                      $<$<COMPILE_LANGUAGE:Fortran>:${SIMINT_Fortran_FLAGS}>)
target_link_libraries(simint PRIVATE ${SIMINT_LINK_FLAGS})

if(${BUILD_SHARED_LIBS})
    target_link_libraries(simint PRIVATE ${LIBC_INTERJECT})

    if(APPLE)
        set_target_properties(simint PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()
endif()

# Installation
# Include directory used for export target
target_include_directories(simint INTERFACE
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS simint
        EXPORT "simintTargets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp")

