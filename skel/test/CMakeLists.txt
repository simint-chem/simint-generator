macro(ADDTESTS TNAME TFILES)
  add_executable(${TNAME} ${TFILES}
                 $<TARGET_OBJECTS:test_valeev>
                 $<TARGET_OBJECTS:test_erd>
                 $<TARGET_OBJECTS:test_libint2>
                 $<TARGET_OBJECTS:test_simint>
                 $<TARGET_OBJECTS:test_common>)

  target_link_libraries(${TNAME} PRIVATE simint
                                 ${LIBERD_LIBRARY} ${LIBINT2_LIBRARY})
                                 
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    target_link_libraries(${TNAME} PRIVATE ifcore)
  else()
    target_link_libraries(${TNAME} PRIVATE gfortran)
  endif()
    

  # Add extra libraries for native MIC
  if(SIMINT_MIC)
    target_link_libraries(${TNAME} PRIVATE intlc imf)
  endif()

  target_include_directories(${TNAME} PRIVATE ${TESTS_CXX_INCLUDES})
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

if(TESTS_BENCHMARK_VALIDATE)
  list(APPEND TESTS_CXX_FLAGS "-DBENCHMARK_VALIDATE")
endif()

if(TESTS_LIBINT2_SIMD)
  list(APPEND TESTS_CXX_FLAGS "-DTESTS_LIBINT2_SIMD")
  list(APPEND TESTS_LIBINT2_FILES "Libint2_simd.cpp")
else()
  list(APPEND TESTS_LIBINT2_FILES "Libint2_scalar.cpp")
endif()


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
  #list(APPEND TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  message(FATAL_ERROR "Unsupported compiler (for now)")
endif()



##################################################
# Libint/liberd comparison
##################################################
set(LIBINT2_PATH "" CACHE FILEPATH "Path to libint2 directory")
set(LIBERD_PATH "" CACHE FILEPATH "Path to liberd lib directory")

if(NOT(LIBINT2_PATH STREQUAL ""))
  list(APPEND TESTS_CXX_INCLUDES "${LIBINT2_PATH}/include/libint2")
  set(LIBINT2_LIBRARY "${LIBINT2_PATH}/lib/libint2.a")
else()
  message(FATAL_ERROR "Error - set libint2 path!")
endif()
if(LIBERD_PATH STREQUAL "")
  message(FATAL_ERROR "Error - set liberd path!")
endif()


##################################################
# OpenMP dependency
##################################################
list(APPEND TESTS_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
list(APPEND TESTS_LINK_FLAGS "${OpenMP_CXX_FLAGS}") 


list(APPEND TESTS_CXX_FLAGS "-DTESTS_USE_LIBERD")
set(LIBERD_LIBRARY "${LIBERD_PATH}/lib/liberd.a")


###########################
# Add OBJECT targets
###########################
# Add common object files
add_library(test_common OBJECT common.cpp)
target_include_directories(test_common PRIVATE ${TESTS_CXX_INCLUDES})
target_compile_options(test_common PRIVATE ${TESTS_CXX_FLAGS})


# Add simint object files
add_library(test_simint OBJECT simint.cpp)
target_include_directories(test_simint PRIVATE ${TESTS_CXX_INCLUDES})
target_compile_options(test_simint PRIVATE ${TESTS_CXX_FLAGS})

# Add valeev object files
add_library(test_valeev OBJECT valeev.cpp)
target_include_directories(test_valeev PRIVATE ${TESTS_CXX_INCLUDES})
target_compile_options(test_valeev PRIVATE ${TESTS_CXX_FLAGS})

# liberd
add_library(test_erd OBJECT ERD.cpp)
target_include_directories(test_erd PRIVATE ${TESTS_CXX_INCLUDES})
target_compile_options(test_erd PRIVATE ${TESTS_CXX_FLAGS})

# Libint2
add_library(test_libint2 OBJECT ${TESTS_LIBINT2_FILES})
target_include_directories(test_libint2 PRIVATE ${TESTS_CXX_INCLUDES})
target_compile_options(test_libint2 PRIVATE ${TESTS_CXX_FLAGS})


ADDTESTS(test_eri test_eri.cpp)

ADDTESTS(benchmark_eri benchmark_eri.cpp)
ADDTESTS(benchmark_erd benchmark_erd.cpp)
ADDTESTS(benchmark_libint2 benchmark_libint2.cpp)
