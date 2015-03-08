add_library(simint_eri STATIC 
                              eri.c
                              eri_ssss.c
           )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
