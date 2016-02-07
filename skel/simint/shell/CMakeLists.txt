add_library(simint_shell OBJECT 
                              shell.c
                              shell_constants.c
           )

target_include_directories(simint_shell PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint_shell PRIVATE ${SIMINT_C_FLAGS})
