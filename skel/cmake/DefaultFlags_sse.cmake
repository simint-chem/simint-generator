if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

  list(APPEND SIMINT_C_FLAGS "-mssse3")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-mssse3")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  list(APPEND SIMINT_C_FLAGS "-mssse3")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-mssse3")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
