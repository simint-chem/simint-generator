if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

  list(APPEND SIMINT_C_FLAGS "-xCORE-AVX512;-fma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-xCORE-AVX512;-fma")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  list(APPEND SIMINT_C_FLAGS "-mavx512f;-mavx512cd;-mfma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-mavx512f;avx512cd;-mfma")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
