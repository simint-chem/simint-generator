if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")

  list(APPEND SIMINT_C_FLAGS "-march=armv8.2-a+sve;-msve-vector-bits=512")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-march=armv8.2-a+sve;-msve-vector-bits=512")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  message(FATAL_ERROR "Unsupported compiler")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()
