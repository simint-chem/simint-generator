add_library(simint_eri STATIC 
                              shell.c
                              shell_constants.c
                              eri_ssss.c
                              eri_ssss_cheby.c
                              eri_ssss_combined.c
                              eri_ssss_taylor.c
                              eri_ssss_taylorcombined.c
                              eri_ssss_split.c
           )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
