file(GLOB SIMINT_FO_FILES FO/*.c)
file(GLOB SIMINT_FOFLAT_FILES FO_flat/*.c)
file(GLOB SIMINT_VREF_FILES vref/*.c)
file(GLOB SIMINT_VREFFLAT_FILES vref_flat/*.c)
#file(GLOB SIMINT_SPLIT_FILES split/*.c)

add_library(simint_eri STATIC 
                              shell.c
                              shell_constants.c
                              vrr/vrr.c
                              hrr/hrr_bra.c
                              hrr/hrr_ket.c
                              ${SIMINT_FO_FILES}
                              ${SIMINT_FOFLAT_FILES}
                              ${SIMINT_VREF_FILES}
                              ${SIMINT_VREFFLAT_FILES}
           )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES "${SIMINT_C_INCLUDES}")
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
