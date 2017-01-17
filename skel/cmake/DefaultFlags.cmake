if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   869 : parameter "xxx" never referenced (ok, but try compiling without it sometimes)
  #   383 : value copied to temporary, reference to temporary used

  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_INTEL")
  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_C_FLAGS "-wd10397;-wd981;-wd869")

  # Temporarily supress unused variable warnings
  # (hope to be fixed in the future)
  list(APPEND SIMINT_C_FLAGS "-wd177")

  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-wd10397;-wd981;-wd383")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-restrict")

  list(APPEND SIMINT_C_LINK_FLAGS "-mkl;-static-intel")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_CLANG")
  elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_GCC")
  endif()

  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-parameter")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-variable")

  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-parameter")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-variable")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-variable")

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()


# SIMINT_VECTOR required to be set (it may be empty here, hence the "")
string(TOLOWER "${SIMINT_VECTOR}" SIMINT_VECTOR_LOWER)

if(SIMINT_VECTOR_LOWER STREQUAL "avx512")
  include(cmake/DefaultFlags_avx512.cmake)
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_AVX512")
elseif(SIMINT_VECTOR_LOWER STREQUAL "micavx512")
  include(cmake/DefaultFlags_micavx512.cmake)
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_MICAVX512")
elseif(SIMINT_VECTOR_LOWER STREQUAL "avxfma")
  include(cmake/DefaultFlags_avxfma.cmake)
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_AVX")
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_FMA")
elseif(SIMINT_VECTOR_LOWER STREQUAL "avx")
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_AVX")
  include(cmake/DefaultFlags_avx.cmake)
elseif(SIMINT_VECTOR_LOWER STREQUAL "sse")
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_SSE")
  include(cmake/DefaultFlags_sse.cmake)
elseif(SIMINT_VECTOR_LOWER STREQUAL "scalar")
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_SCALAR")
  include(cmake/DefaultFlags_scalar.cmake)
else()
  message(FATAL_ERROR "Unsupported vectorization type \"${SIMINT_VECTOR}\"")
endif()

