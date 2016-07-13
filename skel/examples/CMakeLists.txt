###################
# Examples flags
###################
list(APPEND EXAMPLES_C_INCLUDES "${CMAKE_SOURCE_DIR}")

if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  list(APPEND EXAMPLES_C_FLAGS "-std=c99")
else()
  list(APPEND EXAMPLES_C_FLAGS "-Drestrict=__restrict__")
  list(APPEND EXAMPLES_C_FLAGS "-std=c99")
  list(APPEND EXAMPLES_C_FLAGS "-Wall;-Wextra;-pedantic")
  list(APPEND EXAMPLES_LINK_FLAGS m)
endif()



###########################
# Add examples
###########################
add_executable(example1 example1.c)
target_compile_options(example1 PRIVATE ${EXAMPLES_C_FLAGS})
target_include_directories(example1 PRIVATE ${EXAMPLES_C_INCLUDES})

list(APPEND EXAMPLES_LINK_FLAGS simint)
target_link_libraries(example1 PRIVATE ${EXAMPLES_LINK_FLAGS})

