add_executable(compare_f 
                         compare_f.c
              )

target_link_libraries(compare_f simint_boys simint_valeev)

set_target_properties(compare_f PROPERTIES INCLUDE_DIRECTORIES ${SIMINT_C_INCLUDES})
set_target_properties(compare_f PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})
set_target_properties(compare_f PROPERTIES LINK_FLAGS ${SIMINT_LINK_FLAGS})
