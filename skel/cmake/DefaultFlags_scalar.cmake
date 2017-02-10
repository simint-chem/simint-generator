if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

  list(APPEND SIMINT_C_FLAGS "-no-simd;-no-vec")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-no-simd;-no-vec")

endif()
