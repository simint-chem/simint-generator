file(GLOB SIMINT_FO_FILES FO/*.c)
file(GLOB SIMINT_SPLIT_FILES split/*.c)

add_library(simint_eri STATIC 
                              shell.c
                              shell_constants.c
                              eri_split.c
                              eri_splitcombined.c
                              ${SIMINT_FO_FILES}
                              ${SIMINT_SPLIT_FILES}
           )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES "${SIMINT_C_INCLUDES}")
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
