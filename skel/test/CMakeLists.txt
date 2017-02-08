macro(ADDTEST TNAME TFILES)
  add_executable(${TNAME} ${TFILES} $<TARGET_OBJECTS:test_common>)

  target_link_libraries(${TNAME} PRIVATE simint ${SIMINT_TESTS_CXX_LIBRARIES})
                                 
  target_include_directories(${TNAME} PRIVATE ${SIMINT_TESTS_CXX_INCLUDES})
  target_include_directories(${TNAME} SYSTEM PRIVATE ${SIMINT_TESTS_CXX_SYSTEM_INCLUDES})
  target_compile_options(${TNAME} PRIVATE ${SIMINT_TESTS_CXX_FLAGS})

  if( NOT ("${SIMINT_TESTS_LINK_FLAGS}" STREQUAL "") )
    target_link_libraries(${TNAME} PRIVATE ${SIMINT_TESTS_LINK_FLAGS})
  endif()
endmacro() 


###################
# OPTIONS
###################
set(SIMINT_TESTS_BENCHMARK_VALIDATE FALSE CACHE BOOL "Validate integrals in benchmarks")


###################
# Test stuff flags
###################
list(APPEND SIMINT_TESTS_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")
list(APPEND SIMINT_TESTS_CXX_INCLUDES "${CMAKE_BINARY_DIR}/stage")

if(SIMINT_TESTS_BENCHMARK_VALIDATE)
    list(APPEND SIMINT_TESTS_CXX_FLAGS "-DBENCHMARK_VALIDATE;-g")
endif()

message(STATUS "       SIMINT_TESTS_CXX_INCLUDES: ${SIMINT_TESTS_CXX_INCLUDES}")
message(STATUS "SIMINT_TESTS_CXX_SYSTEM_INCLUDES: ${SIMINT_TESTS_CXX_SYSTEM_INCLUDES}")
message(STATUS "          SIMINT_TESTS_CXX_FLAGS: ${SIMINT_TESTS_CXX_FLAGS}")
message(STATUS "      SIMINT_TESTS_CXX_LIBRARIES: ${SIMINT_TESTS_CXX_LIBRARIES}")
message(STATUS "     SIMINT_TESTS_CXX_LINK_FLAGS: ${SIMINT_TESTS_CXX_LINK_FLAGS}")



###########################
# Add OBJECT targets
###########################
# Add common object files
add_library(test_common OBJECT Common.cpp
                               ValeevRef.cpp
                               ${SIMINT_TESTS_CXX_ADDITIONAL_SOURCES})

target_include_directories(test_common PRIVATE ${SIMINT_TESTS_CXX_INCLUDES})
target_include_directories(test_common SYSTEM PRIVATE ${SIMINT_TESTS_CXX_SYSTEM_INCLUDES})
target_compile_options(test_common PRIVATE ${SIMINT_TESTS_CXX_FLAGS})


# Various testing and benchmarking
ADDTEST(test_eri test_eri.cpp)
ADDTEST(test_eri_deriv1 test_eri_deriv1.cpp)
ADDTEST(benchmark_eri benchmark_eri.cpp)
ADDTEST(screentest screentest.cpp)
