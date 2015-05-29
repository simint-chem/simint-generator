macro(ADDTESTS TNAME TFILES)
  add_executable(${TNAME} ${TFILES}
                 $<TARGET_OBJECTS:test_valeev>
                 $<TARGET_OBJECTS:test_erd>
                 $<TARGET_OBJECTS:test_common>)

  target_link_libraries(${TNAME} simint_eri simint_boys
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


###################
# Test stuff flags
###################
list(APPEND TESTS_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")
list(APPEND TESTS_CXX_FLAGS "-std=c++11")
list(APPEND TESTS_CXX_FLAGS "-O3")
list(APPEND TESTS_CXX_FLAGS "-Wall -Wunused-variable -pedantic -Werror -pedantic-errors")
list(APPEND TESTS_CXX_FLAGS "-ggdb")

if(TESTS_BENCHMARK_VALIDATE)
  list(APPEND TESTS_CXX_FLAGS "-DBENCHMARK_VALIDATE")
endif()


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  list(APPEND TESTS_CXX_FLAGS "-qopt-report=4 -wd10397 -restrict")
else()
  list(APPEND TESTS_CXX_FLAGS "-qopt-report=4 -wd10397 -Drestrict=__restrict__")
endif()


##################################################
# AVX, etc
##################################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if(SIMINT_AVX)
    list(APPEND TESTS_CXX_FLAGS "-xavx -DTESTS_AVX -DSIMINT_AVX")
    list(APPEND TESTS_LINK_FLAGS "-xavx")
  elseif(SIMINT_MIC)
    list(APPEND TESTS_CXX_FLAGS "-mmic -DTESTS_MIC -DSIMINT_MIC")
    list(APPEND TESTS_LINK_FLAGS "-mmic")
  endif()
endif()
 

##################################################
# Libint/liberd comparison
##################################################
set(LIBINT2_PATH "" CACHE FILEPATH "Path to libint2 directory")
set(LIBERD_PATH "" CACHE FILEPATH "Path to liberd lib directory")

if(NOT(LIBINT2_PATH STREQUAL ""))
  list(APPEND TESTS_CXX_FLAGS "-DTESTS_USE_LIBINT2")
  list(APPEND TESTS_CXX_INCLUDES "${LIBINT2_PATH}/include/libint2")
  set(LIBINT2_LIBRARY "${LIBINT2_PATH}/lib/libint2.a")
endif()
if(LIBERD_PATH STREQUAL "")
  message(FATAL_ERROR "Error - set liberd path!")
endif()


##################################################
# Optional OpenMP dependency
##################################################
find_package(OpenMP)
if(OPENMP_FOUND)
  message(STATUS "Found OpenMP. Remember to link to OpenMP libraries when building your final project (if needed)")
  list(APPEND TESTS_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
  list(APPEND TESTS_LINK_FLAGS "${OpenMP_CXX_FLAGS}") 
else()
    message(WARNING "OpenMP was not found! Disabling...")
endif()


list(APPEND TESTS_CXX_FLAGS "-DTESTS_USE_LIBERD")
set(LIBERD_LIBRARY "${LIBERD_PATH}/lib/liberd.a")

# Turn lists into a string
# Except for includes. They are handled differently in cmake...
string(REPLACE ";" " " TESTS_CXX_FLAGS "${TESTS_CXX_FLAGS}")
string(REPLACE ";" " " TESTS_LINK_FLAGS "${TESTS_LINK_FLAGS}")

# Add common object files
add_library(test_common OBJECT common.cpp)
set_target_properties(test_common PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_common PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

# Add valeev object files
add_library(test_valeev OBJECT valeev.cpp)
set_target_properties(test_valeev PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_valeev PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})

# Erd interface as well
add_library(test_erd OBJECT erd_interface.cpp)
set_target_properties(test_erd PROPERTIES INCLUDE_DIRECTORIES "${TESTS_CXX_INCLUDES}")
set_target_properties(test_erd PROPERTIES COMPILE_FLAGS ${TESTS_CXX_FLAGS})


ADDTESTS(gen_ref gen_ref.cpp)

ADDTESTS(test_f test_f.cpp)
ADDTESTS(test_eri test_eri.cpp)

ADDTESTS(benchmark_f benchmark_f.cpp)
ADDTESTS(benchmark_FO benchmark_FO.cpp)
ADDTESTS(benchmark_erd benchmark_erd.cpp)
