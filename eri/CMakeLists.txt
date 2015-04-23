add_library(simint_eri STATIC 
                              shell.c
                              shell_constants.c
                              eri_split.c
                              eri_splitcombined.c
                              eri_taylorcombined.c
                              eri_FOcombined.c
           )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES "${SIMINT_C_INCLUDES}")
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
