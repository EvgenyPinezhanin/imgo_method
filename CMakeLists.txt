cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME global_optimization_method)
project(${PROJECT_NAME})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin/release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin/debug)

set(PROJECT_LIBRARY "${PROJECT_NAME}")
set(PROJECT_LIBRARY_PROBLEMS "${PROJECT_NAME}_problems")
set(PROJECT_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(PROJECT_PROBLEMS "${CMAKE_CURRENT_SOURCE_DIR}/GCGen")
set(PROJECT_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src")

include_directories("${PROJECT_INCLUDE}")
include_directories("${PROJECT_PROBLEMS}/sample_src")
include_directories("${PROJECT_PROBLEMS}/src")

# OUTPUT DIR
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set(OUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})
else()
  set(OUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
endif()

# COPY SCRIPT
set(COPY_RESOURCE copy_resource)
add_custom_target(${COPY_RESOURCE} ALL)
add_custom_command(
  TARGET ${COPY_RESOURCE} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_SOURCE_DIR}/scripts/
          ${OUT_DIR}/scripts)

# MAKE DIR
set(MAKE_DIR make_directory)
add_custom_target(${MAKE_DIR} ALL)
add_custom_command(
  TARGET ${MAKE_DIR} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
          ${OUT_DIR}/output_data)

# BUILD
add_subdirectory(GCGen)
add_subdirectory(src)
add_subdirectory(samples)

# REPORT
message( STATUS "")
message( STATUS "General configuration for ${PROJECT_NAME}")
message( STATUS "======================================")
message( STATUS "")
message( STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message( STATUS "")