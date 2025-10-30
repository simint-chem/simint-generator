if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")

  list(APPEND SIMINT_C_FLAGS "-march=armv8.2-a+sve;-msve-vector-bits=512")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-march=armv8.2-a+sve;-msve-vector-bits=512")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Fujitsu")

  message(FATAL_ERROR "Fujitsu is an Unsupported compiler")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  message(FATAL_ERROR "Clang is an Unsupported compiler")

else()

  message(FATAL_ERROR "Unsupported compiler=${CMAKE_C_COMPILER_ID}")

endif()
