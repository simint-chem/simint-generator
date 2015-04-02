file(GLOB LIBERD_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*.f")

add_library(simint_erd 
            ${LIBERD_SOURCE}
            )

##############################################
# Panache library compiler flags, etc
##############################################
# Compiler definitions list need to be split
set_target_properties(simint_erd PROPERTIES COMPILE_FLAGS ${LIBERD_F90_FLAGS})


add_library(simint_erdinterface STATIC 
                              erd_interface.c
           )

set_target_properties(simint_erdinterface PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(simint_erdinterface PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})


# For now, install all header files
#install(TARGETS simint_erd ARCHIVE DESTINATION lib)
