if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #         seems to have lots of false positives, even when it should be obvious
  #   869 : parameter "xxx" never referenced (ok, but try compiling without it sometimes)
  #   383 : value copied to temporary, reference to temporary used

  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_C_FLAGS "-wd10397;-wd981;-wd869")
  list(APPEND SIMINT_C_FLAGS "-DSIMINT_INTEL")

  # Temporarily supress unused variable warnings
  # (hope to be fixed in the future)
  list(APPEND SIMINT_C_FLAGS "-wd177")

  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-wd10397;-wd981;-wd383")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-restrict")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-DSIMINT_INTEL")

  list(APPEND SIMINT_C_LINK_FLAGS "-mkl;-static-intel")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")

  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-parameter")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-variable")
  list(APPEND SIMINT_C_FLAGS "-DSIMINT_GCC")

  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-parameter")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-variable")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-variable")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-DSIMINT_GCC")

else()
  message(FATAL_ERROR "Unsupported compiler")
endif()
