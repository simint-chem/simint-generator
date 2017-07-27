if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

  list(APPEND SIMINT_C_FLAGS "-xAVX")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-xAVX")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  list(APPEND SIMINT_C_FLAGS "-mavx")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-mavx")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
