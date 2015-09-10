file(GLOB SIMINT_ERI_FILES gen/*.c)

add_library(simint_eri OBJECT ${SIMINT_ERI_FILES} )

set_target_properties(simint_eri PROPERTIES INCLUDE_DIRECTORIES "${SIMINT_C_INCLUDES}")
set_target_properties(simint_eri PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
