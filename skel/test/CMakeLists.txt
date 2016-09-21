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

if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   383 : value copied to temporary, reference to temporary used
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-qopt-report=5;-w3;-wd10397;-wd981;-wd383")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-restrict")

else()
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wall;-Wextra;-pedantic")
endif()

if(SIMINT_TESTS_BENCHMARK_VALIDATE)
    list(APPEND SIMINT_TESTS_CXX_FLAGS "-DBENCHMARK_VALIDATE;-g")
endif()


##################################################
# OpenMP dependency
##################################################
list(APPEND SIMINT_TESTS_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
list(APPEND SIMINT_TESTS_LINK_FLAGS "${OpenMP_CXX_FLAGS}") 



##################################################
# Comparison with other packages
##################################################
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
                               Simint.cpp
                               ValeevRef.cpp
                               ${SIMINT_TESTS_CXX_ADDITIONAL_SOURCES})

target_include_directories(test_common PRIVATE ${SIMINT_TESTS_CXX_INCLUDES})
target_include_directories(test_common SYSTEM PRIVATE ${SIMINT_TESTS_CXX_SYSTEM_INCLUDES})
target_compile_options(test_common PRIVATE ${SIMINT_TESTS_CXX_FLAGS})


# Various testing and benchmarking
ADDTEST(test_eri test_eri.cpp)
ADDTEST(benchmark_eri benchmark_eri.cpp)
