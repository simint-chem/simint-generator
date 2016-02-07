add_library(simint_boys OBJECT 
                               boys_init.c
                               boys_longfac.c
                               boys_chebygrid.c
                               boys_shortgrid.c
           )

target_include_directories(simint_boys PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint_boys PRIVATE ${SIMINT_C_FLAGS})
