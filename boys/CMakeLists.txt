add_library(simint_boys STATIC 
                               boys.c
                               boys_shortgrid.c
                               boys_longfac.c
                               boys_cheby.c
           )

set_target_properties(simint_boys PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(simint_boys PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
