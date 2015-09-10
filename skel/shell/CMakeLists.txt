add_library(simint_shell OBJECT 
                              shell.c
                              shell_constants.c
           )

set_target_properties(simint_shell PROPERTIES INCLUDE_DIRECTORIES "${SIMINT_C_INCLUDES}")
set_target_properties(simint_shell PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
