macro(ADDTEST TNAME TFILES)
  add_executable(${TNAME} ${TFILES} $<TARGET_OBJECTS:test_common>)

  target_link_libraries(${TNAME} PRIVATE simint ${TESTS_CXX_LIBRARIES})
                                 
  if(TESTS_ENABLE_LIBERD AND "${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    target_link_libraries(${TNAME} PRIVATE ifcore)
  else()
    target_link_libraries(${TNAME} PRIVATE gfortran m)
  endif()
    

  # Add extra libraries for native MIC
  if(SIMINT_MIC)
    target_link_libraries(${TNAME} PRIVATE intlc imf)
  endif()

  target_include_directories(${TNAME} PRIVATE ${TESTS_CXX_INCLUDES})
  target_include_directories(${TNAME} SYSTEM PRIVATE ${TESTS_CXX_SYSTEM_INCLUDES})
  target_compile_options(${TNAME} PRIVATE ${TESTS_CXX_FLAGS})

  if( NOT ("${TESTS_LINK_FLAGS}" STREQUAL "") )
    target_link_libraries(${TNAME} PRIVATE ${TESTS_LINK_FLAGS})
  endif()
endmacro() 


###################
# OPTIONS
###################
set(TESTS_BENCHMARK_VALIDATE FALSE CACHE BOOL "Validate integrals in benchmarks")
set(TESTS_LIBINT2_SIMD FALSE CACHE BOOL "Libint2 Compiled with SIMD support")


###################
# Test stuff flags
###################
list(APPEND TESTS_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")

if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   383 : value copied to temporary, reference to temporary used
  list(APPEND TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND TESTS_CXX_FLAGS "-qopt-report=5;-w3;-wd10397;-wd981;-wd383")
  list(APPEND TESTS_CXX_FLAGS "-restrict")

else()
  list(APPEND TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  list(APPEND TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND TESTS_CXX_FLAGS "-Wall;-Wextra;-pedantic")
endif()

if(TESTS_BENCHMARK_VALIDATE)
    list(APPEND TESTS_CXX_FLAGS "-DBENCHMARK_VALIDATE;-g")
endif()


##################################################
# OpenMP dependency
##################################################
list(APPEND TESTS_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
list(APPEND TESTS_LINK_FLAGS "${OpenMP_CXX_FLAGS}") 



##################################################
# Comparison with other packages
##################################################
set(LIBINT2_PATH "" CACHE FILEPATH "Path to libint2 directory")
set(LIBERD_PATH "" CACHE FILEPATH "Path to liberd lib directory")

set(TESTS_ENABLE_LIBINT2 FALSE)
set(TESTS_ENABLE_LIBERD FALSE)

# Did the user specify the path to libint2 ?
if(NOT(LIBINT2_PATH STREQUAL ""))
  # Add include path / library
  list(APPEND TESTS_CXX_SYSTEM_INCLUDES "${LIBINT2_PATH}/include/libint2")
  list(APPEND TESTS_CXX_LIBRARIES "${LIBINT2_PATH}/lib/libint2.a")

  # Set some other variables
  set(TESTS_ENABLE_LIBINT2 TRUE)
  list(APPEND TESTS_CXX_FLAGS "-DTESTS_ENABLE_LIBINT2")

  if(TESTS_LIBINT2_SIMD)
    list(APPEND TESTS_CXX_FLAGS "-DTESTS_LIBINT2_SIMD")
    list(APPEND TESTS_CXX_ADDITIONAL_SOURCES "Libint2_simd.cpp")
  else()
    list(APPEND TESTS_CXX_ADDITIONAL_SOURCES "Libint2_scalar.cpp")
  endif()
endif()


if(NOT(LIBERD_PATH STREQUAL ""))
  set(TESTS_ENABLE_LIBERD TRUE)
  list(APPEND TESTS_CXX_FLAGS "-DTESTS_ENABLE_LIBERD")
  list(APPEND TESTS_CXX_ADDITIONAL_SOURCES "ERD.cpp")
  list(APPEND TESTS_CXX_LIBRARIES "${LIBERD_PATH}/lib/liberd.a")
endif()





###########################
# Add OBJECT targets
###########################
# Add common object files
add_library(test_common OBJECT Common.cpp
                               Simint.cpp
                               ValeevRef.cpp
                               ${TESTS_CXX_ADDITIONAL_SOURCES})

target_include_directories(test_common PRIVATE ${TESTS_CXX_INCLUDES})
target_include_directories(test_common SYSTEM PRIVATE ${TESTS_CXX_SYSTEM_INCLUDES})
target_compile_options(test_common PRIVATE ${TESTS_CXX_FLAGS})


# Various testing and benchmarking
ADDTEST(test_eri test_eri.cpp)
ADDTEST(benchmark_eri benchmark_eri.cpp)
#ADDTEST(plotnprim_eri plotnprim_eri.cpp)

if(TESTS_ENABLE_LIBERD)
  ADDTEST(benchmark_erd benchmark_erd.cpp)
  #ADDTEST(plotnprim_erd plotnprim_erd.cpp)
endif()

if(TESTS_ENABLE_LIBINT2)
  ADDTEST(benchmark_libint2 benchmark_libint2.cpp)
endif()
