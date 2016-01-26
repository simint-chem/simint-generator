file(GLOB SIMINT_ERI_FILES gen/*.c)

add_library(simint_eri OBJECT ${SIMINT_ERI_FILES} )

target_include_directories(simint_eri PRIVATE ${SIMINT_C_INCLUDES})
target_compile_options(simint_eri PRIVATE ${SIMINT_C_FLAGS})

