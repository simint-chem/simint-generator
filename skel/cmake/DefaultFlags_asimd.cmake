if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")

  list(APPEND SIMINT_C_FLAGS "-march=native")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-march=native")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  message(FATAL_ERROR "Unsupported compiler")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
