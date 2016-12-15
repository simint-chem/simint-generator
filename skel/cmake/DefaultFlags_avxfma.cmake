include(cmake/DefaultFlags_avx.cmake)

if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

  list(APPEND SIMINT_C_FLAGS "-fma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-fma")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

    list(APPEND SIMINT_C_FLAGS "-mfma")
    list(APPEND SIMINT_TESTS_CXX_FLAGS "-mfma")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
