add_library(simint_eri STATIC 
                              shell.c
                              shell_constants.c
                              eri_erf.c
                              eri_erf_combined.c
                              eri_erf_split.c
                              eri_cheby.c
                              eri_taylor.c
                              eri_taylorcombined.c
           )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
