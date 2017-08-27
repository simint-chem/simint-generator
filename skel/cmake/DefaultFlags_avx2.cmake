if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

  list(APPEND SIMINT_C_FLAGS "-xCORE-AVX2;-fma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-xCORE-AVX2;-fma")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  list(APPEND SIMINT_C_FLAGS "-mavx2;-mfma")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-mavx2;-mfma")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
