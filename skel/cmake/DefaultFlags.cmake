
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  ##########################################################################################
  # Globally-disabled diagnostics:
  # 10397 : Remark about creating .optrpt files
  #   981 : "operands are evaluated in unspecified order"
  #   383 : value copied to temporary, reference to temporary used
  #   869 : parameter "xxx" was never referenced
  ##########################################################################################

  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-qopenmp")
  list(APPEND SIMINT_C_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_C_FLAGS "-wd10397;-wd981;-wd869;-wd177")  # Remove wd177 (unused variable) at some point
  list(APPEND SIMINT_C_LINK_FLAGS "-qopenmp")

  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11;-restrict")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-qopenmp")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-qopt-report=5;-w3")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-wd10397;-wd981;-wd383")
  list(APPEND SIMINT_TESTS_CXX_LINK_FLAGS "-qopenmp")


elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" OR
       "${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

  list(APPEND SIMINT_C_FLAGS "-std=c99")
  list(APPEND SIMINT_C_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-parameter")
  list(APPEND SIMINT_C_FLAGS "-Wno-unused-variable")

  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Drestrict=__restrict__")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-std=c++11")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-parameter")
  list(APPEND SIMINT_TESTS_CXX_FLAGS "-Wno-unused-variable")

  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    # OpenMP is enabled by default with clang
    list(APPEND SIMINT_C_FLAGS "-fopenmp")
    list(APPEND SIMINT_C_LINK_FLAGS "-fopenmp")
    list(APPEND SIMINT_TESTS_CXX_FLAGS "-fopenmp")
    list(APPEND SIMINT_TESTS_CXX_LINK_FLAGS "-fopenmp")
  endif()

else()

  message(FATAL_ERROR "Unsupported compiler")

endif()


# SIMINT_VECTOR required to be set (it may be empty here, hence the "")
string(TOLOWER "${SIMINT_VECTOR}" SIMINT_VECTOR_LOWER)
string(TOUPPER "${SIMINT_VECTOR}" SIMINT_VECTOR_UPPER)

# Handle scalar
if(SIMINT_VECTOR_LOWER STREQUAL "scalar")
  list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_SCALAR")
  include(cmake/DefaultFlags_scalar.cmake)
else()
    if(SIMINT_VECTOR_LOWER MATCHES "scalar-")
        # scalar code, but with other compile flags
        list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_SCALAR")
        string(SUBSTRING ${SIMINT_VECTOR_LOWER} 7 -1 SIMINT_VECTOR_NOSCALAR)
        include(cmake/DefaultFlags_${SIMINT_VECTOR_NOSCALAR}.cmake)
        include(cmake/DefaultFlags_scalar.cmake)
    else()
        list(APPEND SIMINT_CONFIG_DEFINES "SIMINT_${SIMINT_VECTOR_UPPER}")
        include(cmake/DefaultFlags_${SIMINT_VECTOR_LOWER}.cmake)
    endif()
endif()


