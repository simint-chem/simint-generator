add_library(simint_valeev STATIC 
                                 valeev.c
           )

set_target_properties(simint_valeev PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(simint_valeev PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
