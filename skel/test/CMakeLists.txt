macro(ADDTESTS TNAME TFILES)
  add_executable(${TNAME} ${TFILES}
                 $<TARGET_OBJECTS:test_valeev>
                 $<TARGET_OBJECTS:test_erd>
                 $<TARGET_OBJECTS:test_libint2>
                 $<TARGET_OBJECTS:test_libint2_c>
                 $<TARGET_OBJECTS:test_common>)

  target_link_libraries(${TNAME} simint_eri simint_boys simint_shell
                                 ${LIBERD_LIBRARY} ${LIBINT2_LIBRARY})
                                 
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    target_link_libraries(${TNAME} ifcore)
  else()
    target_link_libraries(${TNAME} gfortran)
  endif()
    

  # Add extra libraries for native MIC
  if(SIMINT_MIC)
    target_link_libraries(${TNAME} intlc imf)
  endif()

  set_target_properties(${TNAME} PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
  set_target_properties(${TNAME} PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

  if( NOT ("${TESTS_LINK_FLAGS}" STREQUAL "") )
    set_target_properties(${TNAME} PROPERTIES LINK_FLAGS ${TESTS_LINK_FLAGS})
  endif()
endmacro() 


###################
# OPTIONS
###################
set(TESTS_BENCHMARK_VALIDATE FALSE CACHE BOOL "Validate integrals in benchmarks")
set(PROC_CYCLES_PER_SECOND "" CACHE STRING "Processor cycles per second (in GHz)")


###################
# Test stuff flags
###################
list(APPEND TESTS_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")

if(TESTS_BENCHMARK_VALIDATE)
  list(APPEND TESTS_CXX_FLAGS "-DBENCHMARK_VALIDATE")
endif()


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   383 : value copied to temporary, reference to temporary used
  list(APPEND TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND TESTS_CXX_FLAGS "-qopt-report=5 -w3 -wd10397 -wd981 -wd383")
  list(APPEND TESTS_CXX_FLAGS "-restrict")
  list(APPEND TESTS_CXX_FLAGS "-O3 -ggdb")
else()
  #list(APPEND TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  message(FATAL_ERROR "Unsupported compiler (for now)")
endif()


if("${PROC_CYCLES_PER_SECOND}" STREQUAL "")
  message(WARNING "Warning - set PROC_CYCLES_PER_SECOND for accurate timing in seconds (ticks should still be accurate)")
  set(PROC_CYCLES_PER_SECOND "2.7")
endif()

list(APPEND TESTS_CXX_FLAGS "-DPROC_CYCLES_PER_SECOND=${PROC_CYCLES_PER_SECOND}")



##################################################
# AVX, etc
##################################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if(SIMINT_SSE)
    list(APPEND TESTS_CXX_FLAGS "-msse2 -DSIMINT_SIMD_LEN=2")
    list(APPEND TESTS_LINK_FLAGS "-msse2")
  elseif(SIMINT_AVX)
    list(APPEND TESTS_CXX_FLAGS "-xavx -DSIMINT_SIMD_LEN=4")
    list(APPEND TESTS_LINK_FLAGS "-xavx")
  elseif(SIMINT_MIC)
    list(APPEND TESTS_CXX_FLAGS "-mmic -DSIMINT_SIMD_LEN=8")
    list(APPEND TESTS_LINK_FLAGS "-mmic")
  endif()

  # --no-vec and stuff make no sense here, right?
  if(SIMINT_SCALAR)
    list(APPEND TESTS_CXX_FLAGS "-DSIMINT_SIMD_LEN=1")
  endif()
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

# Turn lists into a string
# Except for includes. They are handled differently in cmake...
string(REPLACE ";" " " TESTS_CXX_FLAGS "${TESTS_CXX_FLAGS}")
string(REPLACE ";" " " TESTS_LINK_FLAGS "${TESTS_LINK_FLAGS}")


###########################
# Add OBJECT targets
###########################
# Add common object files
add_library(test_common OBJECT common.cpp)
set_target_properties(test_common PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_common PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

# Add valeev object files
add_library(test_valeev OBJECT valeev.cpp)
set_target_properties(test_valeev PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_valeev PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

# liberd
add_library(test_erd OBJECT ERD.cpp)
set_target_properties(test_erd PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_erd PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

# Libint2
add_library(test_libint2 OBJECT Libint2.cpp)
set_target_properties(test_libint2 PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_libint2 PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

# Libint2 - c portion
# This is split out so that it can be compiled with the same flags as simint
add_library(test_libint2_c OBJECT libint2/libint2.c)
set_target_properties(test_libint2_c PROPERTIES INCLUDE_DIRECTORIES "${SIMINT_C_INCLUDES};${LIBINT2_PATH}/include/libint2")
set_target_properties(test_libint2_c PROPERTIES COMPILE_FLAGS ${SIMINT_C_FLAGS})


#ADDTESTS(gen_ref gen_ref.cpp)

#ADDTESTS(test_f test_f.cpp)
ADDTESTS(test_eri test_eri.cpp)

#ADDTESTS(benchmark_f benchmark_f.cpp)
ADDTESTS(benchmark_eri benchmark_eri.cpp)
ADDTESTS(benchmark_erd benchmark_erd.cpp)
ADDTESTS(benchmark_libint2 benchmark_libint2.cpp)
