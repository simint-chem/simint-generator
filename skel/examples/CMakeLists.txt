###################
# Examples flags
###################
list(APPEND EXAMPLES_C_INCLUDES "${CMAKE_SOURCE_DIR}")

if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  list(APPEND EXAMPLES_C_FLAGS "-std=c99")
else()
  #list(APPEND EXAMPLES_C_FLAGS "-Drestrict=__restrict__")
  message(FATAL_ERROR "Unsupported compiler (for now)")
endif()



###########################
# Add examples
###########################
add_executable(example1 example1.c)
target_compile_options(example1 PRIVATE ${EXAMPLES_C_FLAGS})
target_include_directories(example1 PRIVATE ${EXAMPLES_C_INCLUDES})
target_link_libraries(example1 PRIVATE simint)

